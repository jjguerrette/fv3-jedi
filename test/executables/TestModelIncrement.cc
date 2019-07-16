/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Utilities/FV3JEDITraits.h"
#include "fv3jedi/Run/RunFV3JEDI.h"
#include "test/base/ModelIncrement.h"

int main(const int argc, const char ** argv) {
  fv3jedi::RunFV3JEDI run(argc, argv);
  test::ModelIncrement<fv3jedi::FV3JEDITraits> tests;
  run.execute(tests);
  return 0;
}

