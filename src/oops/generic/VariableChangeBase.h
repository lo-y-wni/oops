/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_VARIABLECHANGEBASE_H_
#define OOPS_BASE_VARIABLECHANGEBASE_H_

#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// Base class for generic variable transform

template <typename MODEL>
class VariableChangeBase : public util::Printable,
                           private boost::noncopyable {
  typedef Geometry<MODEL>            Geometry_;                           
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  VariableChangeBase() {}
  virtual ~VariableChangeBase() {}

  virtual void linearize(const State_ &, const Geometry_ &) =0;

  virtual void transform(const Increment_ &, Increment_ &) const =0;
  virtual void transformInverse(const Increment_ &, Increment_ &) const =0;
  virtual void transformAdjoint(const Increment_ &, Increment_ &) const =0;
  virtual void transformAdjointInverse(const Increment_ &, Increment_ &) const =0;

 private:
  virtual void print(std::ostream &) const =0;
};

// -----------------------------------------------------------------------------

/// VariableChangeFactory Factory
template <typename MODEL>
class VariableChangeFactory {
  typedef Geometry<MODEL>   Geometry_;
  typedef State<MODEL>      State_;
 public:
  static VariableChangeBase<MODEL> * create(const Geometry_ &, const Variables &,
                                            const eckit::Configuration &, const State_ &);
  virtual ~VariableChangeFactory() { getMakers().clear(); }
 protected:
  explicit VariableChangeFactory(const std::string &);
 private:
  virtual VariableChangeBase<MODEL> * make(const Geometry_ &, const eckit::Configuration &) = 0;
  static std::map < std::string, VariableChangeFactory<MODEL> * > & getMakers() {
    static std::map < std::string, VariableChangeFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class LinearModelMaker : public VariableChangeFactory<MODEL> {
  typedef Geometry<MODEL>   Geometry_;
  virtual VariableChangeBase<MODEL> * make(const Geometry_ & resol, const Variables & vars,
                                           const eckit::Configuration & conf, const State_ & xb)
    { return new T(resol.geometry(), vars, conf, xb.state()); }
 public:
  explicit LinearModelMaker(const std::string & name) : VariableChangeFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
VariableChangeFactory<MODEL>::VariableChangeFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in the variable change factory factory." << std::endl;
    ABORT("Element already registered in VariableChangeFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
VariableChangeBase<MODEL>* VariableChangeFactory<MODEL>::create(const Geometry_ & resol, const Variables & vars,
 77                                            const eckit::Configuration & conf, const State_ & xb) {
  Log::trace() << "VariableChangeBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("version");
  typename std::map<std::string, VariableChangeFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    Log::error() << id << " does not exist in the variable change factory factory." << std::endl;
    ABORT("Element does not exist in VariableChangeFactory.");
  }
  VariableChangeBase<MODEL> * ptr = jerr->second->make(resol, vars, conf, xb);
  Log::trace() << "VariableChangeBase<MODEL>::create done" << std::endl;
  return ptr;
}

// =============================================================================

}  // namespace oops

#endif  // OOPS_BASE_VARIABLECHANGEBASE_H_
