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
printf("Zellen Global  %d %d \n",mypartitioning.nCellsGlobal()[0],mypartitioning.nCellsGlobal()[1]);
printf("Gridsize %d %d \n",mypartitioning.getSubGridSize()[0],mypartitioning.getSubGridSize()[1]);

//settings.nCells=mypartitioning.nCells();

//std::cout << "ncells" << mypartitioning.nCells()[0] << mypartitioning.nCells()[1] << std::endl;

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


// OutputWriterTextParallel myOutputWriterText(myDiscretization, mypartitioning);
OutputWriterParaviewParallel myOutputWriterParaview(myDiscretization, mypartitioning);

int Iterationszahl=0;
// initialize time
double current_time=0;
 //write after initialization

myDiscretization->setBorderVelocityParalell(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
myDiscretization->updateBoundaryFGParalell();
// myOutputWriterParaview.writeFile(current_time);
// myOutputWriterText.writeFile(current_time);
double time_check = 1;
int make_output = 0;

while (current_time<settings.endTime && Iterationszahl< settings.maximumNumberOfIterations )
{
  myDiscretization->updateDeltaT();
  if (current_time+myDiscretization->getDeltaT() > time_check)
  {
    myDiscretization->setDeltaT(current_time+myDiscretization->getDeltaT() - time_check);
    time_check+=1;
    make_output=1;
  }
  
  current_time+=myDiscretization->getDeltaT();
  myDiscretization->calculation();
  myPressureSolver->calculateRHS();
  myPressureSolver->calculateP();
  myDiscretization->updateVelocity();
  //myDiscretization->setBorderVelocity(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
  myDiscretization->setBorderVelocityParalell(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
  myDiscretization->updateBoundaryFGParalell();
  

  if (make_output)
  {
    myOutputWriterParaview.writeFile(current_time);
    make_output=0;
  }
  

  // myOutputWriterText.writeFile(current_time);
Iterationszahl=Iterationszahl+1;
//printf("Iterationszahl %d \n",Iterationszahl);
}
//std::cout<< "Noetige Iterationen " << Iterationszahl <<std::endl;

time1=clock()-tstart;
time1=time1/CLOCKS_PER_SEC;
//printf("Rank %d Iterationszahl %d Laufzeit in s %f\n",mypartitioning.ownRankNo(),Iterationszahl,time1);

//MPI_Finalize();

MPI_Finalize(); // wenn ich () stopt er bei mir das Programm nicht

//std::cout<< "Laufzeit in s " << time1  <<std::endl;
//printf("Laufzeit in s %f \n",time1); //sendet es hier auch mehrmals
}


