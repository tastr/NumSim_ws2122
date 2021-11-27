// !!! Alte main   !!!!

#include "computation.h"


Computation::Computation(Settings settings)
:settings_(settings)
{

}


void Computation::runSimulation(Settings settings)
{
// testing 
// MyTestFunctions myTest;
// assert(myTest.testFunction2(settings));
double time1, tstart;
tstart=clock();

//create objects of classes
MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);



std::shared_ptr<Discretization> myDiscretization;
Partitioning mypartitioning(settings);
settings.nCells=mypartitioning.nCells();

std::cout << "ncells" << settings.nCells[0] << settings.nCells[1] << std::endl;
if (settings.useDonorCell)
{
  myDiscretization=std::make_shared<DonorCell>(settings,mypartitioning);
} else
{
  myDiscretization=std::make_shared<CentralDifferences>(settings,mypartitioning);
}
std::shared_ptr<PressureSolver> myPressureSolver;
if (settings.pressureSolver == "SOR")
{
  myPressureSolver= std::make_shared<SOR>(*myDiscretization); //reference is used
} else
{
  myPressureSolver= std::make_shared<GaussSeidel>(*myDiscretization); //reference is used
}


// std::shared_ptr<Discretization> pointer_to_myDiscretization (& myDiscretization); //vermutlich gibt es da einen besseren Weg, aber den habe ich nicht gefunden...
OutputWriterText myOutputWriterText(myDiscretization, mypartitioning);
OutputWriterParaview myOutputWriterParaview(myDiscretization, mypartitioning);
int Iterationszahl=0;
// initialize time
double current_time=0;
 //write after initialization

myDiscretization->setBorderVelocity(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
myDiscretization->updateBoundaryFG();
myOutputWriterParaview.writeFile(current_time);
myOutputWriterText.writeFile(current_time);

while (current_time<settings.endTime && Iterationszahl< settings.maximumNumberOfIterations )
{
  myDiscretization->updateDeltaT();
  current_time+=myDiscretization->getDeltaT();
  myDiscretization->calculation();
  myPressureSolver->calculateRHS();
  myPressureSolver->calculateP();
  myDiscretization->updateVelocity();
  myDiscretization->setBorderVelocity(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
  myDiscretization->updateBoundaryFG();
  
  myOutputWriterParaview.writeFile(current_time);
  myOutputWriterText.writeFile(current_time);
Iterationszahl=Iterationszahl+1;
}
std::cout<< "Noetige Iterationen " << Iterationszahl <<std::endl;
time1=clock()-tstart;
time1=time1/CLOCKS_PER_SEC;
MPI_Finalize();
std::cout<< "Laufzeit in s " << time1  <<std::endl;
}

