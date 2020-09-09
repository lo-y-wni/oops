/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_GETKFSOLVER_H_
#define OOPS_ASSIMILATION_GETKFSOLVER_H_

#include <Eigen/Dense>
#include <cfloat>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/gletkfInterface.h"
#include "oops/assimilation/LETKFSolverParameters.h"
#include "oops/assimilation/LocalEnsembleSolver.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/StateEnsemble.h"
#include "oops/generic/VerticalLocEV.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/util/Logger.h"

namespace oops {

/*!
 * An implementation of the GETKF from Lei 2018 JAMES
 *
 * Lei, L., Whitaker, J. S., & Bishop, C. ( 2018). Improving assimilation 
 * of radiance observations by implementing model space localization in an 
 * ensemble Kalman filter. Journal of Advances in Modeling Earth Systems, 10, 
 * 3221– 3232. https://doi.org/10.1029/2018MS001468
 */
template <typename MODEL, typename OBS>
class GETKFSolver : public LocalEnsembleSolver<MODEL, OBS> {
  typedef Departures<OBS>           Departures_;
  typedef DeparturesEnsemble<OBS>   DeparturesEnsemble_;
  typedef Geometry<MODEL>           Geometry_;
  typedef GeometryIterator<MODEL>   GeometryIterator_;
  typedef IncrementEnsemble<MODEL>  IncrementEnsemble_;
  typedef ObsEnsemble<OBS>          ObsEnsemble_;
  typedef ObsErrors<OBS>            ObsErrors_;
  typedef Observations<OBS>         Observations_;
  typedef ObsSpaces<OBS>            ObsSpaces_;
  typedef State4D<MODEL>            State4D_;
  typedef StateEnsemble<MODEL>      StateEnsemble_;
  typedef VerticalLocEV<MODEL>      VerticalLocEV_;

 public:
  /// Constructor (allocates Wa, wa, HZb_,
  /// saves options from the config, computes VerticalLocEV_)
  GETKFSolver(ObsSpaces_ &, const Geometry_ &, const eckit::Configuration &, size_t);

  Observations_ computeHofX(const StateEnsemble_ &, size_t) override;

  /// entire KF update (computeWeights+applyWeights) for a grid point GeometryIterator_
  void measurementUpdate(const IncrementEnsemble_ &, const GeometryIterator_ &,
                         IncrementEnsemble_ &) override;

 private:
  /// Computes weights
  void computeWeights(const Departures_ &, const DeparturesEnsemble_ &,
                      const DeparturesEnsemble_ &, const ObsErrors_ &);

  /// Applies weights and adds posterior inflation
  void applyWeights(const IncrementEnsemble_ &, IncrementEnsemble_ &,
                    const GeometryIterator_ &);

 private:
  LETKFSolverParameters options_;

  // parameters
  size_t nens_;
  const Geometry_ & geometry_;
  VerticalLocEV_ vertloc_;
  size_t neig_;
  size_t nanal_;

  DeparturesEnsemble_ HZb_;

  Eigen::MatrixXd Wa_;  // transformation matrix for ens. perts. Xa_=Xf*Wa
  Eigen::VectorXd wa_;  // transformation matrix for ens. mean xa_=xf*wa
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
GETKFSolver<MODEL, OBS>::GETKFSolver(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                const eckit::Configuration & config, size_t nens)
  : LocalEnsembleSolver<MODEL, OBS>(obspaces, geometry, config, nens),
    nens_(nens), geometry_(geometry),
    vertloc_(geometry_, config.getSubConfiguration("letkf.vertical localization")),
    neig_(vertloc_.neig()), nanal_(neig_*nens_), HZb_(obspaces, nanal_)
{
  options_.deserialize(config);
  const LETKFInflationParameters & inflopt = options_.infl;

  // parse inflation options
  Log::info() << "Multiplicative inflation multCoeff=" <<
                 inflopt.mult << std::endl;

  if (inflopt.dortpp()) {
      Log::info() << "RTPP inflation will be applied with rtppCoeff=" <<
                    inflopt.rtpp << std::endl;
  } else {
      Log::info() << "RTPP inflation is not applied rtppCoeff is out of bounds (0,1], rtppCoeff="
                  << inflopt.rtpp << std::endl;
  }

  if (inflopt.dortps()) {
    Log::info() << "RTPS inflation will be applied with rtpsCoeff=" <<
                    inflopt.rtps << std::endl;
  } else {
    Log::info() << "RTPS inflation is not applied rtpsCoeff is out of bounds (0,1], rtpsCoeff="
                << inflopt.rtps << std::endl;
  }

  // pre-allocate transformation matrices
  Wa_.resize(nanal_, nens);
  wa_.resize(nanal_);
}

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
Observations<OBS> GETKFSolver<MODEL, OBS>::computeHofX(const StateEnsemble_ & ens_xx,
                                                       size_t iteration) {
  // compute H(x) for the original ensemble members
  Observations_ yb_mean = LocalEnsembleSolver<MODEL, OBS>::computeHofX(ens_xx, iteration);

  // modulate ensemble of obs
  State4D_ xx_mean(ens_xx.mean());
  IncrementEnsemble_ dx(ens_xx, xx_mean, xx_mean[0].variables());
  IncrementEnsemble_ Ztmp(geometry_, xx_mean[0].variables(), ens_xx[0].validTimes(), neig_);
  size_t ii = 0;
  for (size_t iens = 0; iens < nens_; ++iens) {
    vertloc_.modulateIncrement(dx[iens], Ztmp);
    for (size_t ieig = 0; ieig < neig_; ++ieig) {
      State4D_ tmpState = xx_mean;
      tmpState += Ztmp[ieig];
      Observations_ tmpObs = this->hofx_.compute(tmpState);
      HZb_[ii++] = tmpObs - yb_mean;
    }
  }

  return yb_mean;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GETKFSolver<MODEL, OBS>::computeWeights(const Departures_ & dy,
                                             const DeparturesEnsemble_ & Yb,
                                             const DeparturesEnsemble_ & YbOrig,
                                             const ObsErrors_ & R) {
  // compute transformation matrix, save in Wa_, wa_
  // Yb(nobs,neig*nens), YbOrig(nobs,nens)
  // uses GSI GETKF code

  const int nobsl = dy.nobs();

  // cast oops objects to eigen<double> obejects
  // then cast eigen<doble> to eigen<float>
  Eigen::MatrixXd edy = dy.packEigen();
  Eigen::MatrixXf edy_f = edy.cast<float>();

  Eigen::MatrixXd eYb = Yb.packEigen();
  Eigen::MatrixXf eYb_f = eYb.cast<float>();

  Eigen::MatrixXd eYb2 = YbOrig.packEigen();
  Eigen::MatrixXf eYb2_f = eYb2.cast<float>();

  Eigen::MatrixXd eR = R.packInverseVarianceEigen();
  Eigen::MatrixXf eR_f = eR.cast<float>();

  Eigen::MatrixXf Wa_f(this->nanal_, this->nens_);
  Eigen::VectorXf wa_f(this->nanal_);

  // call into GSI interface to compute Wa and wa
  const int getkf_inflation = 0;
  const int denkf = 0;
  const int getkf = 1;
  letkf_core_f90(nobsl, eYb_f.data(), eYb2_f.data(), edy_f.data(),
                 wa_f.data(), Wa_f.data(),
                 eR_f.data(), nanal_, neig_,
                 getkf_inflation, denkf, getkf);
  this->Wa_ = Wa_f.cast<double>();
  this->wa_ = wa_f.cast<double>();
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GETKFSolver<MODEL, OBS>::applyWeights(const IncrementEnsemble_ & bkg_pert,
                                           IncrementEnsemble_ & ana_pert,
                                           const GeometryIterator_ & i) {
  // apply Wa_, wa_

  const LETKFInflationParameters & inflopt = options_.infl;

  // allocate tmp arrays
  LocalIncrement gptmpl = bkg_pert[0][0].getLocal(i);
  std::vector<double> tmp = gptmpl.getVals();
  size_t ngp = tmp.size();
  Eigen::MatrixXd Xb(ngp, nens_);

  // loop through analysis times and ens. members
  for (unsigned itime=0; itime < bkg_pert[0].size(); ++itime) {
    // cast bkg_pert ensemble at grid point i as an Eigen matrix Xb
    // modulates Xb
    Xb = vertloc_.modulateIncrement(bkg_pert, i, itime);

    // postmulptiply
    Eigen::VectorXd xa = Xb*wa_;  // ensemble mean update
    Eigen::MatrixXd Xa = Xb*Wa_;  // ensemble perturbation update

    // compute non-modulated Xb for RTPP and RTPS
    if (inflopt.dortpp() || inflopt.dortps()) {
      Xb.resize(ngp, nens_);
      for (size_t iens=0; iens < nens_; ++iens) {
        LocalIncrement gp = bkg_pert[iens][itime].getLocal(i);
        std::vector<double> tmp1 = gp.getVals();
        for (size_t iv=0; iv < ngp; ++iv) {
          Xb(iv, iens) = tmp1[iv];
        }
      }
    }

    // RTPP inflation
    if (inflopt.dortpp()) {
      Xa = (1-inflopt.rtpp)*Xa+inflopt.rtpp*Xb;
    }

    // RTPS inflation
    double eps = DBL_EPSILON;
    if (inflopt.dortps()) {
      // posterior spread
      Eigen::ArrayXd asprd = Xa.array().square().rowwise().sum()/(nens_-1);
      asprd = asprd.sqrt();
      asprd = (asprd < eps).select(eps, asprd);  // avoid nan overflow for vars with no spread

      // prior spread
      Eigen::ArrayXd fsprd = Xb.array().square().rowwise().sum()/(nens_-1);
      fsprd = fsprd.sqrt();
      fsprd = (fsprd < eps).select(eps, fsprd);

      // rtps inflation factor
      Eigen::ArrayXd rtpsInfl = inflopt.rtps*((fsprd-asprd)/asprd) + 1;
      rtpsInfl = (rtpsInfl < inflopt.rtpsInflMin()).select(inflopt.rtpsInflMin(), rtpsInfl);
      rtpsInfl = (rtpsInfl > inflopt.rtpsInflMax()).select(inflopt.rtpsInflMax(), rtpsInfl);

      // inlfate perturbation matrix
      Xa.array().colwise() *= rtpsInfl;
    }

    // assign Xa_ to ana_pert
    for (size_t iens=0; iens < nens_; ++iens) {
      for (size_t iv=0; iv < ngp; ++iv) {
        tmp[iv] = Xa(iv, iens)+xa(iv);
      }
      gptmpl.setVals(tmp);
      ana_pert[iens][itime].setLocal(gptmpl, i);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GETKFSolver<MODEL, OBS>::measurementUpdate(const IncrementEnsemble_ & bkg_pert,
                                                const GeometryIterator_ & i,
                                                IncrementEnsemble_ & ana_pert) {
  // create the local subset of observations
  ObsSpaces_ local_obs(this->obspaces_, *i, this->obsconf_);
  Departures_ local_omb(local_obs, this->omb_);

  if (local_omb.nobs() == 0) {
    // no obs. so no need to update Wa_ and wa_
    // ana_pert[i]=bkg_pert[i]
    this->copyLocalIncrement(bkg_pert, i, ana_pert);
  } else {
    // if obs are present do normal KF update
    DeparturesEnsemble_ local_Yb(local_obs, this->Yb_);
    DeparturesEnsemble_ local_HZ(local_obs, HZb_);
    // create local obs errors
    ObsErrors_ local_R(this->obsconf_, local_obs);
    computeWeights(local_omb, local_HZ, local_Yb, local_R);
    applyWeights(bkg_pert, ana_pert, i);
  }
}

}  // namespace oops
#endif  // OOPS_ASSIMILATION_GETKFSOLVER_H_