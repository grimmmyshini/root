#include "RooFit/TestStatistics/MinuitFcnAD.h"
#include "LikelihoodSerial.h"

#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooAbsPdf.h"
#include "RooNaNPacker.h"

#include <fstream>
#include <iomanip>  // std::setprecision

using namespace std;

namespace RooFit {
namespace TestStatistics {

RooArgSet getParameters(RooAbsReal const& funct) {
   RooArgSet out;
   funct.getParameters(nullptr, out);
   return out;
}

MinuitFcnAD::MinuitFcnAD(RooAbsReal* funct, RooMinimizer *context, bool verbose)
   : RooAbsMinimizerFcn(getParameters(*funct), context, verbose), funct_(funct)
{ }

void MinuitFcnAD::Gradient(const double *x, double *grad) const
{
  DoEval(x);
  funct_->evaluateGradient(grad);
}

void MinuitFcnAD::GradientWithPrevResult(const double *x, double *grad, double *previous_grad, double *previous_g2,
                                           double *previous_gstep) const
{
  Gradient(x, grad);
}

double MinuitFcnAD::DoEval(const double *x) const {

  // Set the parameter values for this iteration
  for (unsigned index = 0; index < _nDim; index++) {
    if (_logfile) (*_logfile) << x[index] << " " ;
    SetPdfParamVal(index,x[index]);
  }

  // Calculate the function for these parameters
  RooAbsReal::setHideOffset(kFALSE) ;
  double fvalue = funct_->getVal();
  RooAbsReal::setHideOffset(kTRUE) ;

  if (!std::isfinite(fvalue) || RooAbsReal::numEvalErrors() > 0 || fvalue > 1e30) {
    printEvalErrors();
    RooAbsReal::clearEvalErrorLog() ;
    _numBadNLL++ ;

    if (_doEvalErrorWall) {
      const double badness = RooNaNPacker::unpackNaN(fvalue);
      fvalue = (std::isfinite(_maxFCN) ? _maxFCN : 0.) + _recoverFromNaNStrength * badness;
    }
  } else {
    if (_evalCounter > 0 && _evalCounter == _numBadNLL) {
      // This is the first time we get a valid function value; while before, the
      // function was always invalid. For invalid  cases, we returned values > 0.
      // Now, we offset valid values such that they are < 0.
      _funcOffset = -fvalue;
    }
    fvalue += _funcOffset;
    _maxFCN = std::max(fvalue, _maxFCN);
  }

  // Optional logging
  if (_logfile)
    (*_logfile) << setprecision(15) << fvalue << setprecision(4) << endl;
  if (_verbose) {
    cout << "\nprevFCN" << (funct_->isOffsetting()?"-offset":"") << " = " << setprecision(10)
         << fvalue << setprecision(4) << "  " ;
    cout.flush() ;
  }

  _evalCounter++ ;

  return fvalue;
}

double MinuitFcnAD::DoDerivative(const double * /*x*/, unsigned int /*icoord*/) const
{
   throw std::runtime_error("MinuitFcnAD::DoDerivative is not implemented, please use Gradient instead.");
}

} // namespace TestStatistics
} // namespace RooFit
