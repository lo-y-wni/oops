/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsBias.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "eckit/config/Configuration.h"
#include "model/ObsBiasIncrement.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
ObsBias::ObsBias(const ObsSpaceQG &, const eckit::Configuration & conf)
  : bias_(ntypes, 0.0), active_(false), geovars_(), hdiags_() {
  oops::Log::info() << "ObsBias: conf = " << conf << std::endl;
  eckit::LocalConfiguration biasconf;
  if (conf.has("ObsBias")) {
    conf.get("ObsBias", biasconf);
    active_ = biasconf.has("stream") || biasconf.has("uwind") ||
              biasconf.has("vwind") || biasconf.has("wspeed");
  }
  if (active_) {
    if (biasconf.has("stream")) bias_[0] = biasconf.getDouble("stream");
    if (biasconf.has("uwind"))  bias_[1] = biasconf.getDouble("uwind");
    if (biasconf.has("vwind"))  bias_[2] = biasconf.getDouble("vwind");
    if (biasconf.has("wspeed")) bias_[3] = biasconf.getDouble("wspeed");
    std::string strn = "";
    for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
      if (jj > 0) strn += ", ";
      std::ostringstream strs;
      strs << bias_[jj];
      strn += strs.str();
    }
    oops::Log::info() << "ObsBias::ObsBias created, bias = " << strn << std::endl;
  }
}
// -----------------------------------------------------------------------------
ObsBias::ObsBias(const ObsBias & other, const bool copy)
  : bias_(ntypes, 0.0), active_(other.active_),
    geovars_(other.geovars_), hdiags_(other.hdiags_)
{
  if (active_ && copy) {
    for (unsigned int jj = 0; jj < ntypes; ++jj) bias_[jj] = other.bias_[jj];
  }
}
// -----------------------------------------------------------------------------
ObsBias & ObsBias::operator+=(const ObsBiasIncrement & dx) {
  if (active_) {
    for (unsigned int jj = 0; jj < ntypes; ++jj) bias_[jj] += dx[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
ObsBias & ObsBias::operator=(const ObsBias & rhs) {
  if (active_) {
    for (unsigned int jj = 0; jj < ntypes; ++jj) bias_[jj] = rhs.bias_[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
double ObsBias::norm() const {
  double zz = 0.0;
  if (active_) {
    double ztmp = 0.0;
    std::size_t ii = 0;
    for (unsigned int jj = 0; jj < ntypes; ++jj) {
      ztmp = bias_[jj]*bias_[jj];
      zz += ztmp;
      if (ztmp > 0.0) ++ii;
    }
    zz = std::sqrt(zz/ii);
  }
  return zz;
}
// -----------------------------------------------------------------------------
void ObsBias::print(std::ostream & os) const {
  if (active_) {
    std::string strn = "";
    for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
      if (jj > 0) strn += ", ";
      std::ostringstream strs;
      strs << bias_[jj];
      strn += strs.str();
    }
    os << std::endl << "ObsBias = " << strn;
  }
}
// -----------------------------------------------------------------------------
}  // namespace qg

