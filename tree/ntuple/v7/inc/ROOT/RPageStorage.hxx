/// \file ROOT/RPageStorage.hxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-07-19
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RPageStorage
#define ROOT7_RPageStorage

#include <ROOT/RCluster.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleMetrics.hxx>
#include <ROOT/RNTupleOptions.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RPage.hxx>
#include <ROOT/RPageAllocator.hxx>
#include <ROOT/RSpan.hxx>
#include <ROOT/RStringView.hxx>

#include <atomic>
#include <cstddef>
#include <deque>
#include <functional>
#include <memory>
#include <shared_mutex>
#include <unordered_set>
#include <vector>

namespace ROOT {
namespace Experimental {

class RNTupleModel;
// TODO(jblomer): factory methods to create tree sinks and sources outside Detail namespace

namespace Detail {

class RColumn;
class RColumnElementBase;
class RNTupleCompressor;
class RNTupleDecompressor;
class RPagePool;
class RFieldBase;

enum class EPageStorageType {
   kSink,
   kSource,
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RPageStorage
\ingroup NTuple
\brief Common functionality of an ntuple storage for both reading and writing

The RPageStore provides access to a storage container that keeps the bits of pages and clusters comprising
an ntuple.  Concrete implementations can use a TFile, a raw file, an object store, and so on.
*/
// clang-format on
class RPageStorage {
public:
   /// The interface of a task scheduler to schedule page (de)compression tasks
   class RTaskScheduler {
   public:
      virtual ~RTaskScheduler() = default;
      /// Start a new set of tasks
      virtual void Reset() = 0;
      /// Take a callable that represents a task
      virtual void AddTask(const std::function<void(void)> &taskFunc) = 0;
      /// Blocks until all scheduled tasks finished
      virtual void Wait() = 0;
   };

   /// A sealed page contains the bytes of a page as written to storage (packed & compressed).  It is used
   /// as an input to UnsealPages() as well as to transfer pages between different storage media.
   /// RSealedPage does _not_ own the buffer it is pointing to in order to not interfere with the memory management
   /// of concrete page sink and page source implementations.
   struct RSealedPage {
      const void *fBuffer = nullptr;
      std::uint32_t fSize = 0;
      std::uint32_t fNElements = 0;

      RSealedPage() = default;
      RSealedPage(const void *b, std::uint32_t s, std::uint32_t n) : fBuffer(b), fSize(s), fNElements(n) {}
      RSealedPage(const RSealedPage &other) = delete;
      RSealedPage& operator =(const RSealedPage &other) = delete;
      RSealedPage(RSealedPage &&other) = default;
      RSealedPage& operator =(RSealedPage &&other) = default;
   };

   using SealedPageSequence_t = std::deque<RSealedPage>;
   /// A range of sealed pages referring to the same column that can be used for vector commit
   struct RSealedPageGroup {
      DescriptorId_t fPhysicalColumnId;
      SealedPageSequence_t::const_iterator fFirst;
      SealedPageSequence_t::const_iterator fLast;

      RSealedPageGroup(DescriptorId_t d, SealedPageSequence_t::const_iterator b, SealedPageSequence_t::const_iterator e)
         : fPhysicalColumnId(d), fFirst(b), fLast(e)
      {
      }
   };

protected:
   std::string fNTupleName;
   RTaskScheduler *fTaskScheduler = nullptr;

public:
   explicit RPageStorage(std::string_view name);
   RPageStorage(const RPageStorage &other) = delete;
   RPageStorage& operator =(const RPageStorage &other) = delete;
   RPageStorage(RPageStorage &&other) = default;
   RPageStorage& operator =(RPageStorage &&other) = default;
   virtual ~RPageStorage();

   /// Whether the concrete implementation is a sink or a source
   virtual EPageStorageType GetType() = 0;

   struct RColumnHandle {
      DescriptorId_t fPhysicalId = kInvalidDescriptorId;
      const RColumn *fColumn = nullptr;

      /// Returns true for a valid column handle; fColumn and fPhysicalId should always either both
      /// be valid or both be invalid.
      explicit operator bool() const { return fPhysicalId != kInvalidDescriptorId && fColumn; }
   };
   /// The column handle identifies a column with the current open page storage
   using ColumnHandle_t = RColumnHandle;

   /// Register a new column.  When reading, the column must exist in the ntuple on disk corresponding to the meta-data.
   /// When writing, every column can only be attached once.
   virtual ColumnHandle_t AddColumn(DescriptorId_t fieldId, const RColumn &column) = 0;
   /// Unregisters a column.  A page source decreases the reference counter for the corresponding active column.
   /// For a page sink, dropping columns is currently a no-op.
   virtual void DropColumn(ColumnHandle_t columnHandle) = 0;

   /// Every page store needs to be able to free pages it handed out.  But Sinks and sources have different means
   /// of allocating pages.
   virtual void ReleasePage(RPage &page) = 0;

   /// Page storage implementations have their own metrics. The RPageSink and RPageSource classes provide
   /// a default set of metrics.
   virtual RNTupleMetrics &GetMetrics() = 0;
   /// Returns the NTuple name.
   const std::string &GetNTupleName() const { return fNTupleName; }

   void SetTaskScheduler(RTaskScheduler *taskScheduler) { fTaskScheduler = taskScheduler; }
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RPageSink
\ingroup NTuple
\brief Abstract interface to write data into an ntuple

The page sink takes the list of columns and afterwards a series of page commits and cluster commits.
The user is responsible to commit clusters at a consistent point, i.e. when all pages corresponding to data
up to the given entry number are committed.
*/
// clang-format on
class RPageSink : public RPageStorage {
private:
   /// Used to map the IDs of the descriptor to the physical IDs issued during header/footer serialization
   Internal::RNTupleSerializer::RContext fSerializationContext;

protected:
   /// Default I/O performance counters that get registered in fMetrics
   struct RCounters {
      RNTupleAtomicCounter &fNPageCommitted;
      RNTupleAtomicCounter &fSzWritePayload;
      RNTupleAtomicCounter &fSzZip;
      RNTupleAtomicCounter &fTimeWallWrite;
      RNTupleAtomicCounter &fTimeWallZip;
      RNTupleTickCounter<RNTupleAtomicCounter> &fTimeCpuWrite;
      RNTupleTickCounter<RNTupleAtomicCounter> &fTimeCpuZip;
   };
   std::unique_ptr<RCounters> fCounters;
   RNTupleMetrics fMetrics;

   std::unique_ptr<RNTupleWriteOptions> fOptions;

   /// Helper to zip pages and header/footer; includes a 16MB (kMAXZIPBUF) zip buffer.
   /// There could be concrete page sinks that don't need a compressor.  Therefore, and in order to stay consistent
   /// with the page source, we leave it up to the derived class whether or not the compressor gets constructed.
   std::unique_ptr<RNTupleCompressor> fCompressor;

   /// Remembers the starting cluster id for the next cluster group
   std::uint64_t fNextClusterInGroup = 0;
   /// Used to calculate the number of entries in the current cluster
   NTupleSize_t fPrevClusterNEntries = 0;
   /// Keeps track of the number of elements in the currently open cluster. Indexed by column id.
   std::vector<RClusterDescriptor::RColumnRange> fOpenColumnRanges;
   /// Keeps track of the written pages in the currently open cluster. Indexed by column id.
   std::vector<RClusterDescriptor::RPageRange> fOpenPageRanges;
   RNTupleDescriptorBuilder fDescriptorBuilder;

   virtual void CreateImpl(const RNTupleModel &model, unsigned char *serializedHeader, std::uint32_t length) = 0;
   virtual RNTupleLocator CommitPageImpl(ColumnHandle_t columnHandle, const RPage &page) = 0;
   virtual RNTupleLocator
   CommitSealedPageImpl(DescriptorId_t physicalColumnId, const RPageStorage::RSealedPage &sealedPage) = 0;
   /// Vector commit of preprocessed pages. The `ranges` array specifies a range of sealed pages to be
   /// committed for each column.  The returned vector contains, in order, the RNTupleLocator for each
   /// page on each range in `ranges`, i.e. the first N entries refer to the N pages in `ranges[0]`,
   /// followed by M entries that refer to the M pages in `ranges[1]`, etc.
   /// The default is to call `CommitSealedPageImpl` for each page; derived classes may provide an
   /// optimized implementation though.
   virtual std::vector<RNTupleLocator> CommitSealedPageVImpl(std::span<RPageStorage::RSealedPageGroup> ranges);
   /// Returns the number of bytes written to storage (excluding metadata)
   virtual std::uint64_t CommitClusterImpl(NTupleSize_t nEntries) = 0;
   /// Returns the locator of the page list envelope of the given buffer that contains the serialized page list.
   /// Typically, the implementation takes care of compressing and writing the provided buffer.
   virtual RNTupleLocator CommitClusterGroupImpl(unsigned char *serializedPageList, std::uint32_t length) = 0;
   virtual void CommitDatasetImpl(unsigned char *serializedFooter, std::uint32_t length) = 0;

   /// Helper for streaming a page. This is commonly used in derived, concrete page sinks. Note that if
   /// compressionSetting is 0 (uncompressed) and the page is mappable, the returned sealed page will
   /// point directly to the input page buffer.  Otherwise, the sealed page references an internal buffer
   /// of fCompressor.  Thus, the buffer pointed to by the RSealedPage should never be freed.
   /// Usage of this method requires construction of fCompressor.
   RSealedPage SealPage(const RPage &page, const RColumnElementBase &element, int compressionSetting);

   /// Seal a page using the provided buffer.
   static RSealedPage SealPage(const RPage &page, const RColumnElementBase &element,
      int compressionSetting, void *buf);

   /// Enables the default set of metrics provided by RPageSink. `prefix` will be used as the prefix for
   /// the counters registered in the internal RNTupleMetrics object.
   /// This set of counters can be extended by a subclass by calling `fMetrics.MakeCounter<...>()`.
   ///
   /// A subclass using the default set of metrics is always responsible for updating the counters
   /// appropriately, e.g. `fCounters->fNPageCommited.Inc()`
   ///
   /// Alternatively, a subclass might provide its own RNTupleMetrics object by overriding the
   /// GetMetrics() member function.
   void EnableDefaultMetrics(const std::string &prefix);

public:
   RPageSink(std::string_view ntupleName, const RNTupleWriteOptions &options);

   RPageSink(const RPageSink&) = delete;
   RPageSink& operator=(const RPageSink&) = delete;
   RPageSink(RPageSink&&) = default;
   RPageSink& operator=(RPageSink&&) = default;
   ~RPageSink() override;

   /// Guess the concrete derived page source from the file name (location)
   static std::unique_ptr<RPageSink> Create(std::string_view ntupleName, std::string_view location,
                                            const RNTupleWriteOptions &options = RNTupleWriteOptions());
   EPageStorageType GetType() final { return EPageStorageType::kSink; }
   /// Returns the sink's write options.
   const RNTupleWriteOptions &GetWriteOptions() const { return *fOptions; }

   ColumnHandle_t AddColumn(DescriptorId_t fieldId, const RColumn &column) final;
   void DropColumn(ColumnHandle_t /*columnHandle*/) final {}

   /// Physically creates the storage container to hold the ntuple (e.g., a keys a TFile or an S3 bucket)
   /// To do so, Create() calls CreateImpl() after updating the descriptor.
   /// Create() associates column handles to the columns referenced by the model
   void Create(RNTupleModel &model);
   /// Write a page to the storage. The column must have been added before.
   void CommitPage(ColumnHandle_t columnHandle, const RPage &page);
   /// Write a preprocessed page to storage. The column must have been added before.
   void CommitSealedPage(DescriptorId_t physicalColumnId, const RPageStorage::RSealedPage &sealedPage);
   /// Write a vector of preprocessed pages to storage. The corresponding columns must have been added before.
   void CommitSealedPageV(std::span<RPageStorage::RSealedPageGroup> ranges);
   /// Finalize the current cluster and create a new one for the following data.
   /// Returns the number of bytes written to storage (excluding meta-data).
   std::uint64_t CommitCluster(NTupleSize_t nEntries);
   /// Write out the page locations (page list envelope) for all the committed clusters since the last call of
   /// CommitClusterGroup (or the beginning of writing).
   void CommitClusterGroup();
   /// Finalize the current cluster and the entrire data set.
   void CommitDataset();

   /// Get a new, empty page for the given column that can be filled with up to nElements.  If nElements is zero,
   /// the page sink picks an appropriate size.
   virtual RPage ReservePage(ColumnHandle_t columnHandle, std::size_t nElements) = 0;

   /// Returns the default metrics object.  Subclasses might alternatively provide their own metrics object by overriding this.
   RNTupleMetrics &GetMetrics() override { return fMetrics; };
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RPageSource
\ingroup NTuple
\brief Abstract interface to read data from an ntuple

The page source is initialized with the columns of interest. Alias columns from projected fields are mapped to the
corresponding physical columns. Pages from the columns of interest can then be mapped into memory.
The page source also gives access to the ntuple's meta-data.
*/
// clang-format on
class RPageSource : public RPageStorage {
public:
   /// An RAII wrapper used for the read-only access to RPageSource::fDescriptor. See GetExclDescriptorGuard().
   class RSharedDescriptorGuard {
      const RNTupleDescriptor &fDescriptor;
      std::shared_mutex &fLock;

   public:
      RSharedDescriptorGuard(const RNTupleDescriptor &desc, std::shared_mutex &lock) : fDescriptor(desc), fLock(lock)
      {
         fLock.lock_shared();
      }
      RSharedDescriptorGuard(const RSharedDescriptorGuard &) = delete;
      RSharedDescriptorGuard &operator=(const RSharedDescriptorGuard &) = delete;
      RSharedDescriptorGuard(RSharedDescriptorGuard &&) = delete;
      RSharedDescriptorGuard &operator=(RSharedDescriptorGuard &&) = delete;
      ~RSharedDescriptorGuard() { fLock.unlock_shared(); }
      const RNTupleDescriptor *operator->() const { return &fDescriptor; }
      const RNTupleDescriptor &GetRef() const { return fDescriptor; }
   };

   /// An RAII wrapper used for the writable access to RPageSource::fDescriptor. See GetSharedDescriptorGuard().
   class RExclDescriptorGuard {
      RNTupleDescriptor &fDescriptor;
      std::shared_mutex &fLock;

   public:
      RExclDescriptorGuard(RNTupleDescriptor &desc, std::shared_mutex &lock) : fDescriptor(desc), fLock(lock)
      {
         fLock.lock();
      }
      RExclDescriptorGuard(const RExclDescriptorGuard &) = delete;
      RExclDescriptorGuard &operator=(const RExclDescriptorGuard &) = delete;
      RExclDescriptorGuard(RExclDescriptorGuard &&) = delete;
      RExclDescriptorGuard &operator=(RExclDescriptorGuard &&) = delete;
      ~RExclDescriptorGuard()
      {
         fDescriptor.IncGeneration();
         fLock.unlock();
      }
      RNTupleDescriptor *operator->() const { return &fDescriptor; }
      void MoveIn(RNTupleDescriptor &&desc) { fDescriptor = std::move(desc); }
   };

private:
   RNTupleDescriptor fDescriptor;
   mutable std::shared_mutex fDescriptorLock;

protected:
   /// Default I/O performance counters that get registered in fMetrics
   struct RCounters {
      RNTupleAtomicCounter &fNReadV;
      RNTupleAtomicCounter &fNRead;
      RNTupleAtomicCounter &fSzReadPayload;
      RNTupleAtomicCounter &fSzReadOverhead;
      RNTupleAtomicCounter &fSzUnzip;
      RNTupleAtomicCounter &fNClusterLoaded;
      RNTupleAtomicCounter &fNPageLoaded;
      RNTupleAtomicCounter &fNPagePopulated;
      RNTupleAtomicCounter &fTimeWallRead;
      RNTupleAtomicCounter &fTimeWallUnzip;
      RNTupleTickCounter<RNTupleAtomicCounter> &fTimeCpuRead;
      RNTupleTickCounter<RNTupleAtomicCounter> &fTimeCpuUnzip;
      RNTupleCalcPerf &fBandwidthReadUncompressed;
      RNTupleCalcPerf &fBandwidthReadCompressed;
      RNTupleCalcPerf &fBandwidthUnzip;
      RNTupleCalcPerf &fFractionReadOverhead;
      RNTupleCalcPerf &fCompressionRatio;
   };

   /// Keeps track of the requested physical column IDs. When using alias columns (projected fields), physical
   /// columns may be requested multiple times.
   class RActivePhysicalColumns {
   private:
      std::vector<DescriptorId_t> fIDs;
      std::vector<std::size_t> fRefCounters;

   public:
      void Insert(DescriptorId_t physicalColumnID);
      void Erase(DescriptorId_t physicalColumnID);
      RCluster::ColumnSet_t ToColumnSet() const;
   };

   std::unique_ptr<RCounters> fCounters;
   /// Wraps the I/O counters and is observed by the RNTupleReader metrics
   RNTupleMetrics fMetrics;

   RNTupleReadOptions fOptions;
   /// The active columns are implicitly defined by the model fields or views
   RActivePhysicalColumns fActivePhysicalColumns;

   /// Helper to unzip pages and header/footer; comprises a 16MB (kMAXZIPBUF) unzip buffer.
   /// Not all page sources need a decompressor (e.g. virtual ones for chains and friends don't), thus we
   /// leave it up to the derived class whether or not the decompressor gets constructed.
   std::unique_ptr<RNTupleDecompressor> fDecompressor;

   virtual RNTupleDescriptor AttachImpl() = 0;
   // Only called if a task scheduler is set. No-op be default.
   virtual void UnzipClusterImpl(RCluster * /* cluster */)
      { }

   /// Helper for unstreaming a page. This is commonly used in derived, concrete page sources.  The implementation
   /// currently always makes a memory copy, even if the sealed page is uncompressed and in the final memory layout.
   /// The optimization of directly mapping pages is left to the concrete page source implementations.
   /// Usage of this method requires construction of fDecompressor.
   std::unique_ptr<unsigned char []> UnsealPage(const RSealedPage &sealedPage, const RColumnElementBase &element);

   /// Enables the default set of metrics provided by RPageSource. `prefix` will be used as the prefix for
   /// the counters registered in the internal RNTupleMetrics object.
   /// A subclass using the default set of metrics is responsible for updating the counters
   /// appropriately, e.g. `fCounters->fNRead.Inc()`
   /// Alternatively, a subclass might provide its own RNTupleMetrics object by overriding the
   /// GetMetrics() member function.
   void EnableDefaultMetrics(const std::string &prefix);

   /// Note that the underlying lock is not recursive. See GetSharedDescriptorGuard() for further information.
   RExclDescriptorGuard GetExclDescriptorGuard() { return RExclDescriptorGuard(fDescriptor, fDescriptorLock); }

public:
   RPageSource(std::string_view ntupleName, const RNTupleReadOptions &fOptions);
   RPageSource(const RPageSource&) = delete;
   RPageSource& operator=(const RPageSource&) = delete;
   RPageSource(RPageSource &&) = delete;
   RPageSource &operator=(RPageSource &&) = delete;
   ~RPageSource() override;
   /// Guess the concrete derived page source from the file name (location)
   static std::unique_ptr<RPageSource> Create(std::string_view ntupleName, std::string_view location,
                                              const RNTupleReadOptions &options = RNTupleReadOptions());
   /// Open the same storage multiple time, e.g. for reading in multiple threads
   virtual std::unique_ptr<RPageSource> Clone() const = 0;

   EPageStorageType GetType() final { return EPageStorageType::kSource; }
   const RNTupleReadOptions &GetReadOptions() const { return fOptions; }

   /// Takes the read lock for the descriptor. Multiple threads can take the lock concurrently.
   /// The underlying std::shared_mutex, however, is neither read nor write recursive:
   /// within one thread, only one lock (shared or exclusive) must be acquired at the same time. This requires special
   /// care in sections protected by GetSharedDescriptorGuard() and GetExclDescriptorGuard() especially to avoid that
   /// the locks are acquired indirectly (e.g. by a call to GetNEntries()).
   /// As a general guideline, no other method of the page source should be called (directly or indirectly) in a
   /// guarded section.
   const RSharedDescriptorGuard GetSharedDescriptorGuard() const
   {
      return RSharedDescriptorGuard(fDescriptor, fDescriptorLock);
   }

   ColumnHandle_t AddColumn(DescriptorId_t fieldId, const RColumn &column) override;
   void DropColumn(ColumnHandle_t columnHandle) override;

   /// Open the physical storage container for the tree
   void Attach() { GetExclDescriptorGuard().MoveIn(AttachImpl()); }
   NTupleSize_t GetNEntries();
   NTupleSize_t GetNElements(ColumnHandle_t columnHandle);
   ColumnId_t GetColumnId(ColumnHandle_t columnHandle);

   /// Allocates and fills a page that contains the index-th element
   virtual RPage PopulatePage(ColumnHandle_t columnHandle, NTupleSize_t globalIndex) = 0;
   /// Another version of PopulatePage that allows to specify cluster-relative indexes
   virtual RPage PopulatePage(ColumnHandle_t columnHandle, const RClusterIndex &clusterIndex) = 0;

   /// Read the packed and compressed bytes of a page into the memory buffer provided by selaedPage. The sealed page
   /// can be used subsequently in a call to RPageSink::CommitSealedPage.
   /// The fSize and fNElements member of the sealedPage parameters are always set. If sealedPage.fBuffer is nullptr,
   /// no data will be copied but the returned size information can be used by the caller to allocate a large enough
   /// buffer and call LoadSealedPage again.
   virtual void
   LoadSealedPage(DescriptorId_t physicalColumnId, const RClusterIndex &clusterIndex, RSealedPage &sealedPage) = 0;

   /// Populates all the pages of the given cluster ids and columns; it is possible that some columns do not
   /// contain any pages.  The page source may load more columns than the minimal necessary set from `columns`.
   /// To indicate which columns have been loaded, LoadClusters() must mark them with SetColumnAvailable().
   /// That includes the ones from the `columns` that don't have pages; otherwise subsequent requests
   /// for the cluster would assume an incomplete cluster and trigger loading again.
   /// LoadClusters() is typically called from the I/O thread of a cluster pool, i.e. the method runs
   /// concurrently to other methods of the page source.
   virtual std::vector<std::unique_ptr<RCluster>> LoadClusters(std::span<RCluster::RKey> clusterKeys) = 0;

   /// Parallel decompression and unpacking of the pages in the given cluster. The unzipped pages are supposed
   /// to be preloaded in a page pool attached to the source. The method is triggered by the cluster pool's
   /// unzip thread. It is an optional optimization, the method can safely do nothing. In particular, the
   /// actual implementation will only run if a task scheduler is set. In practice, a task scheduler is set
   /// if implicit multi-threading is turned on.
   void UnzipCluster(RCluster *cluster);

   /// Returns the default metrics object.  Subclasses might alternatively override the method and provide their own metrics object.
   RNTupleMetrics &GetMetrics() override { return fMetrics; };
};

} // namespace Detail

} // namespace Experimental
} // namespace ROOT

#endif
