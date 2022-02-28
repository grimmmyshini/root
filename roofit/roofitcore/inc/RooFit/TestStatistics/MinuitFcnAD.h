#ifndef ROOT_ROOFIT_TESTSTATISTICS_MinuitFcnAD
#define ROOT_ROOFIT_TESTSTATISTICS_MinuitFcnAD

#include "RooArgList.h"
#include "RooRealVar.h"
#include <RooFit/TestStatistics/RooAbsL.h>
#include <RooFit/TestStatistics/LikelihoodWrapper.h>
#include <RooFit/TestStatistics/LikelihoodGradientWrapper.h>
#include "RooAbsMinimizerFcn.h"

#include <Fit/ParameterSettings.h>
#include <Fit/Fitter.h>
#include "Math/IFunction.h" // ROOT::Math::IMultiGradFunction

// forward declaration
class RooAbsReal;
class RooMinimizer;

namespace RooFit {
namespace TestStatistics {

class MinuitFcnAD : public ROOT::Math::IMultiGradFunction, public RooAbsMinimizerFcn {
public:

   MinuitFcnAD(RooAbsReal* func, RooMinimizer *context, bool verbose = false);

   inline ROOT::Math::IMultiGradFunction *Clone() const override { return new MinuitFcnAD(*this); }

   // used inside Minuit:
   inline bool returnsInMinuit2ParameterSpace() const override { return false /*gradient->usesMinuitInternalValues()*/; }

   inline void setOptimizeConstOnFunction(RooAbsArg::ConstOpCode opcode, Bool_t doAlsoTrackingOpt) override
   { }

private:
   /// IMultiGradFunction override necessary for Minuit
   double DoEval(const double *x) const override;

public:
   /// IMultiGradFunction overrides necessary for Minuit
   void Gradient(const double *x, double *grad) const override;
   void GradientWithPrevResult(const double *x, double *grad, double *previous_grad, double *previous_g2,
                               double *previous_gstep) const override;

   /// Part of IMultiGradFunction interface, used widely both in Minuit and in RooFit.
   inline unsigned int NDim() const override { return _nDim; }

   inline std::string getFunctionName() const override { return funct_->GetName(); }

   inline std::string getFunctionTitle() const override { return funct_->GetTitle(); }

   inline void setOffsetting(Bool_t flag) override { funct_->enableOffsetting(flag); }

private:
   /// This override should not be used in this class, so it throws.
   double DoDerivative(const double *x, unsigned int icoord) const override;

private:
   RooAbsReal *funct_;

 public:
   mutable bool minuit_internal_roofit_x_mismatch_ = false;
};

} // namespace TestStatistics
} // namespace RooFit

#endif // ROOT_ROOFIT_TESTSTATISTICS_MinuitFcnAD
