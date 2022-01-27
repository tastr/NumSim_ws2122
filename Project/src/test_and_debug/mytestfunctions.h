#pragma once

#include "settings.h"
#include "output_writer/output_writer_text.h"
#include "output_writer/output_writer_paraview.h"
#include "pressure_solver/pressuresolver.h"
#include "pressure_solver/sor.h"
#include "pressure_solver/gaussseidel.h"
#include "discretization_storage/discretization.h"
#include "discretization_storage/donorcell.h"
#include "discretization_storage/centraldifferences.h"
#include "test_and_debug/mytestfunctions.h"


class MyTestFunctions
{
    private:

    public:
    bool testFunction1(Settings settings);
    bool testFunction2(Settings settings);

};