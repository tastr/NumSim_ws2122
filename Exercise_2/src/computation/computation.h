#pragma once

#include <iostream>
#include <cstdlib>
#include <cassert>

#include <time.h>
#include <memory>

#include "settings.h"
// #include "output_writer/output_writer_text.h"
// #include "output_writer/output_writer_paraview.h"
#include "../output_writer/output_writer_text_parallel.h"
#include "../output_writer/output_writer_paraview_parallel.h"
#include "pressure_solver/pressuresolver.h"
#include "pressure_solver/sor.h"
#include "pressure_solver/gaussseidel.h"
#include "discretization_storage/discretization.h"
#include "discretization_storage/donorcell.h"
#include "discretization_storage/centraldifferences.h"
#include "partitioning/partitioning.h"
// #include "test_and_debug/mytestfunctions.h"

class Computation
{
protected:
    Settings settings_;

public:
    Computation(Settings settings);
    void runSimulation(Settings settings);
};