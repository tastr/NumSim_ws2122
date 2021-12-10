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
printf("Rank: %d. Zellen Global:  %d %d. Zellen lokal: %d %d. \n", mypartitioning.ownRankNo(),mypartitioning.nCellsGlobal()[0],mypartitioning.nCellsGlobal()[1],mypartitioning.nCells()[0],mypartitioning.nCells()[1]);
// printf("Gridsize %d %d \n",mypartitioning.getSubGridSize()[0],mypartitioning.getSubGridSize()[1]);

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

 // profiler
 double time_computation = 0;
 double time_computation_start = 0;
 double time_communication = 0;
 double delta_communication_pressure = 0;
 double time_communication_pressure = 0;
 double time_communication_start = 0;


myDiscretization->setBorderVelocityParalell(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
myDiscretization->updateBoundaryFGParalell();
// myOutputWriterParaview.writeFile(current_time);
// myOutputWriterText.writeFile(current_time);
double time_check = 1;
double time_stretch = 1.2;
int make_output = 0;

while (current_time<settings.endTime)
{

  myDiscretization->updateDeltaT();
  

   if (current_time+((myDiscretization->getDeltaT())*time_stretch) > time_check)
  {
    myDiscretization->setDeltaT(time_check - current_time);
    time_check+=1;
    make_output=1;
  }
  
  current_time+=myDiscretization->getDeltaT();
  // start stopping the time
  time_computation_start=MPI_Wtime();
  myDiscretization->calculation();
  myPressureSolver->calculateRHS();
  delta_communication_pressure = myPressureSolver->calculateP();
  myDiscretization->updateVelocity();

  time_computation = time_computation+MPI_Wtime()-time_computation_start - delta_communication_pressure;

  time_communication_start = MPI_Wtime();
  //myDiscretization->setBorderVelocity(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
  myDiscretization->setBorderVelocityParalell(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
  myDiscretization->updateBoundaryFGParalell();
  time_communication = time_communication+MPI_Wtime()-time_communication_start;
  time_communication_pressure += delta_communication_pressure;
  

  if (make_output)
  {
    myOutputWriterParaview.writeFile(current_time);
    // myOutputWriterText.writeFile(current_time);

    make_output=0;
  }
  

Iterationszahl=Iterationszahl+1;
//printf("Iterationszahl %d \n",Iterationszahl);
}
//std::cout<< "Noetige Iterationen " << Iterationszahl <<std::endl;

time1=clock()-tstart;
time1=time1/CLOCKS_PER_SEC;
printf("Rank: %d. Iterationszahl: %d. Laufzeit in s %f. Zeit für Berechnung: %f. Zeit für Kommunikation (ohne Druck): %f. Zeit für Druck kommuikation: %f.\n",mypartitioning.ownRankNo(),Iterationszahl,time1, time_computation,time_communication, time_communication_pressure);

MPI_Allreduce(&time_computation,&time_computation,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
MPI_Allreduce(&time_communication,&time_communication,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
MPI_Allreduce(&time_communication_pressure,&time_communication_pressure,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
MPI_Allreduce(&time1,&time1,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

if (mypartitioning.ownRankNo()==0)
{
  printf("!!!Maximum aller Ranks: Laufzeit: %f Zeit für Berechnung: %f. Zeit für Kommunikation (ohne Druck): %f. Zeit für Druck kommuikation: %f.\n", time1, time_computation,time_communication, time_communication_pressure);
} 



//MPI_Finalize();

MPI_Finalize(); // wenn ich () stopt er bei mir das Programm nicht

//std::cout<< "Laufzeit in s " << time1  <<std::endl;
//printf("Laufzeit in s %f \n",time1); //sendet es hier auch mehrmals
}


