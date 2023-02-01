/*
 * Project: RooFit
 * Authors:
 *   Garima Singh, CERN 2022
 *
 * Copyright (c) 2022, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#ifndef RooFit_RooFuncWrapper_h
#define RooFit_RooFuncWrapper_h

#include "TROOT.h"
#include "TSystem.h"
#include "RooAbsReal.h"
#include "RooGlobalFunc.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooListProxy.h"

#include <memory>
#include <string>

/// @brief  A wrapper class to store a C++ function of type 'double (*)(double* )'.
/// The parameters can be accessed as x[<relative position of param in paramSet>] in the function body.
/// @tparam Func Function pointer to the generated function.
template <typename Func = double (*)(double *), typename Grad = void (*)(double *, double *)>
class RooFuncWrapper final : public RooAbsReal {
public:
   RooFuncWrapper(const char *name, const char *title, std::string const &funcBody, RooArgSet const &paramSet)
      : RooAbsReal{name, title}, _params{"!params", "List of parameters", this}
   {
      std::string funcName = name;
      std::string gradName = funcName + "_grad";
      std::string requestName = funcName + "_req";
      std::string wrapperName = funcName + "_derivativeWrapper";

      gInterpreter->Declare("#pragma cling optimize(2)");

      // Declare the function
      std::string bodyWithSig = "double " + funcName + "(double* x) {" + funcBody + "}";
      bool comp = gInterpreter->Declare(bodyWithSig.c_str());
      if (!comp) {
         std::stringstream errorMsg;
         errorMsg << "Function " << funcName << " could not be compiled. See above for details.";
         coutE(InputArguments) << errorMsg.str() << std::endl;
         throw std::runtime_error(errorMsg.str().c_str());
      }
      _func = (Func)gInterpreter->ProcessLine((funcName + ";").c_str());

      // Calculate gradient
      gInterpreter->ProcessLine("#include <Math/CladDerivator.h>");
      // disable clang-format for making the following code unreadable.
      // clang-format off
      std::string requestFunc = "#pragma clad ON\n"
                                 "void " + requestName + "() {\n"
                                 "  clad::gradient(" + funcName + ");\n"
                                 "}\n"
                                 "#pragma clad OFF";
      // clang-format on
      comp = gInterpreter->Declare(requestFunc.c_str());
      if (!comp) {
         std::stringstream errorMsg;
         errorMsg << "Function " << funcName << " could not be be differentiated. See above for details.";
         coutE(InputArguments) << errorMsg.str() << std::endl;
         throw std::runtime_error(errorMsg.str().c_str());
      }

      // Extract parameters
      for (auto *param : paramSet) {
         if (!dynamic_cast<RooAbsReal *>(param)) {
            std::stringstream errorMsg;
            errorMsg << "In creation of function " << funcName
                     << " wrapper: input param expected to be of type RooAbsReal.";
            coutE(InputArguments) << errorMsg.str() << std::endl;
            throw std::runtime_error(errorMsg.str().c_str());
         }
         _params.add(*param);
      }

      // Build a wrapper over the derivative to hide clad specific types such as 'array_ref'.
      // disable clang-format for making the following code unreadable.
      // clang-format off
      std::string nParams = std::to_string(_params.size());
      std::string dWrapper = "void " + wrapperName + "(double* in, double* out) {\n"
                             "  clad::array_ref<double> cladOut(out, " + nParams + ");\n"
                             "  " + gradName + "(in, cladOut);\n"
                             "}";
      // clang-format on
      gInterpreter->Declare(dWrapper.c_str());
      _grad = (Grad)gInterpreter->ProcessLine((wrapperName + ";").c_str());
   }

   RooFuncWrapper(const RooFuncWrapper &other, const char *name = nullptr)
      : RooAbsReal(other, name), _params("!params", this, other._params), _func(other._func)
   {
   }

   TObject *clone(const char *newname) const override { return new RooFuncWrapper(*this, newname); }

   double defaultErrorLevel() const override { return 0.5; }

   void getGradient(double *out) const
   {
      std::vector<double> paramVals;
      getVariablesAsDouble(paramVals);

      _grad(paramVals.data(), out);
   }

protected:
   void getVariablesAsDouble(std::vector<double> &out) const
   {
      out.reserve(_params.size());
      std::transform(_params.begin(), _params.end(), std::back_inserter(out),
                     [](RooAbsArg *obj) { return static_cast<RooAbsReal *>(obj)->getValV(); });
   }

   double evaluate() const override
   {
      std::vector<double> paramVals;
      getVariablesAsDouble(paramVals);

      return _func(paramVals.data());
   }

private:
   RooListProxy _params;
   Func _func;
   Grad _grad;
};

#endif
