/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTFCT3DVAR_H_
#define OOPS_ASSIMILATION_COSTFCT3DVAR_H_

#include "eckit/config/Configuration.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb3D.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// 3D-Var Cost Function
/*!
 * This class is not really necessary since it is only a special
 * case of the more general 4D-Var cost function. It is provided
 * for readability.
 */

// -----------------------------------------------------------------------------

template<typename MODEL> class CostFct3DVar : public CostFunction<MODEL> {
  typedef Increment<MODEL>           Increment_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef ControlVariable<MODEL>     CtrlVar_;
  typedef CostFunction<MODEL>        CostFct_;
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef Model<MODEL>               Model_;

 public:
  CostFct3DVar(const eckit::Configuration &, const Geometry_ &, const Model_ &);
  virtual ~CostFct3DVar() {}

  void runTLM(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>, const bool) const override;
  void runADJ(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>, const bool) const override;
  void zeroAD(CtrlInc_ &) const override;

  void runNL(CtrlVar_ &, PostProcessor<State_>&) const override;

 private:
  void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_>&) const override;

  CostJb3D<MODEL>     * newJb(const eckit::Configuration &, const Geometry_ &,
                              const CtrlVar_ &) const override;
  CostJo<MODEL>       * newJo(const eckit::Configuration &) const override;
  CostTermBase<MODEL> * newJc(const eckit::Configuration &, const Geometry_ &) const override;

  util::Duration windowLength_;
  util::DateTime windowBegin_;
  util::DateTime windowEnd_;
  util::DateTime windowHalf_;
  util::Duration zero_;
  const Variables ctlvars_;
};

// =============================================================================

template<typename MODEL>
CostFct3DVar<MODEL>::CostFct3DVar(const eckit::Configuration & config,
                                  const Geometry_ & resol, const Model_ & model)
  : CostFunction<MODEL>::CostFunction(config, resol, model),
    windowLength_(), windowHalf_(), zero_(0), ctlvars_(config)
{
  Log::trace() << "CostFct3DVar::CostFct3DVar start" << std::endl;
  windowLength_ = util::Duration(config.getString("window_length"));
  windowBegin_ = util::DateTime(config.getString("window_begin"));
  windowEnd_ = windowBegin_ + windowLength_;
  windowHalf_ = windowBegin_ + windowLength_/2;
  this->setupTerms(config);
  Log::trace() << "CostFct3DVar::CostFct3DVar done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostJb3D<MODEL> * CostFct3DVar<MODEL>::newJb(const eckit::Configuration & jbConf,
                                             const Geometry_ & resol,
                                             const CtrlVar_ & xb) const {
  Log::trace() << "CostFct3DVar::newJb" << std::endl;
  ASSERT(xb.state().checkStatesNumber(1));
  return new CostJb3D<MODEL>(jbConf, resol, ctlvars_, zero_, xb.state()[0]);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostJo<MODEL> * CostFct3DVar<MODEL>::newJo(const eckit::Configuration & joConf) const {
  Log::trace() << "CostFct3DVar::newJo" << std::endl;
  return new CostJo<MODEL>(joConf, windowBegin_, windowEnd_, windowLength_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostTermBase<MODEL> * CostFct3DVar<MODEL>::newJc(const eckit::Configuration & jcConf,
                                                 const Geometry_ &) const {
  Log::trace() << "CostFct3DVar::newJc" << std::endl;
// For now there is no Jc that can work with 3D-Var
  CostTermBase<MODEL> * pjc = 0;
  return pjc;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct3DVar<MODEL>::runNL(CtrlVar_ & xx,
                                PostProcessor<State_> & post) const {
  Log::trace() << "CostFct3DVar::runNL start" << std::endl;
  ASSERT(xx.state().checkStatesNumber(1));
  ASSERT(xx.state()[0].validTime() == windowHalf_);
  State_ xm(xx.state()[0].geometry(), CostFct_::getModel().variables(), windowHalf_);
  xm = xx.state()[0];
  CostFct_::getModel().forecast(xm, xx.modVar(), util::Duration(0), post);
  Log::trace() << "CostFct3DVar::runNL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct3DVar<MODEL>::runTLM(CtrlInc_ & dx,
                                 PostProcessorTLAD<MODEL> & cost,
                                 PostProcessor<Increment_> post,
                                 const bool) const {
  Log::trace() << "CostFct3DVar::runTLM start" << std::endl;
  ASSERT(dx.state()[0].validTime() == windowHalf_);
  CostFct_::getTLM().forecastTL(dx.state()[0], dx.modVar(), util::Duration(0), post, cost);
  ASSERT(dx.state()[0].validTime() == windowHalf_);
  Log::trace() << "CostFct3DVar::runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct3DVar<MODEL>::zeroAD(CtrlInc_ & dx) const {
  Log::trace() << "CostFct3DVar::zeroAD start" << std::endl;
  dx.state()[0].zero(windowHalf_);
  dx.modVar().zero();
  dx.obsVar().zero();
  Log::trace() << "CostFct3DVar::zeroAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct3DVar<MODEL>::runADJ(CtrlInc_ & dx,
                                 PostProcessorTLAD<MODEL> & cost,
                                 PostProcessor<Increment_> post,
                                 const bool) const {
  Log::trace() << "CostFct3DVar::runADJ start" << std::endl;
  ASSERT(dx.state()[0].validTime() == windowHalf_);
  CostFct_::getTLM().forecastAD(dx.state()[0], dx.modVar(), util::Duration(0), post, cost);
  ASSERT(dx.state()[0].validTime() == windowHalf_);
  Log::trace() << "CostFct3DVar::runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostFct3DVar<MODEL>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                  PostProcessor<Increment_> &) const {
  Log::trace() << "CostFct3DVar::addIncr start" << std::endl;
  ASSERT(xx.state().checkStatesNumber(1));
  ASSERT(xx.state()[0].validTime() == windowHalf_);
  ASSERT(dx.state()[0].validTime() == windowHalf_);
  xx.state()[0] += dx.state()[0];
  Log::trace() << "CostFct3DVar::addIncr done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCT3DVAR_H_
