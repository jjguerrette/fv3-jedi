/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_GETVALUESTRAJFV3JEDI_H_
#define SRC_GETVALUESTRAJFV3JEDI_H_

#include <ostream>

#include "oops/util/Printable.h"
#include "Fortran.h"

namespace fv3jedi {

// -----------------------------------------------------------------------------

class GetValuesTrajFV3JEDI : public util::Printable {
 public:
  GetValuesTrajFV3JEDI();
  ~GetValuesTrajFV3JEDI();

  int & toFortran() {return keyGetValuesTraj_;}
  const int & toFortran() const {return keyGetValuesTraj_;}

 private:
  void print(std::ostream &) const {}
  F90ootrj keyGetValuesTraj_;
};

// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // SRC_GETVALUESTRAJFV3JEDI_H_
