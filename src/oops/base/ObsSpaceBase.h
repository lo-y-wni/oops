/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSSPACEBASE_H_
#define OOPS_BASE_OBSSPACEBASE_H_

#include <boost/noncopyable.hpp>

#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// Base class for observation spaces.

class ObsSpaceBase : public util::Printable,
                     private boost::noncopyable {
 public:
  ObsSpaceBase(const eckit::Configuration &, const util::DateTime & bgn, const util::DateTime & end)
  : winbgn_(bgn), winend_(end) {}

  virtual ~ObsSpaceBase() {}

/// Access information
  const util::DateTime & windowStart() const {return winbgn_;}
  const util::DateTime & windowEnd() const {return winend_;}

/// Pure virtual methods
  virtual void generateDistribution(const eckit::Configuration &) = 0;

 private:
  const util::DateTime winbgn_;
  const util::DateTime winend_;
};

}  // namespace oops

#endif  // OOPS_BASE_OBSSPACEBASE_H_