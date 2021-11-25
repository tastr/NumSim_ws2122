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
std::shared_ptr<Discretization> myDiscretization;
if (settings.useDonorCell)
{
  myDiscretization=std::make_shared<DonorCell>(settings);
} else
{
  myDiscretization=std::make_shared<CentralDifferences>(settings);
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
// OutputWriterText myOutputWriterText(myDiscretization);
// OutputWriterParaview myOutputWriterParaview(myDiscretization);
int Iterationszahl=0;
// initialize time
double current_time=0;
 //write after initialization

myDiscretization->setBorderVelocity(settings.dirichletBcTop, settings.dirichletBcLeft, settings.dirichletBcRight, settings.dirichletBcBottom);
myDiscretization->updateBoundaryFG();
// myOutputWriterParaview.writeFile(current_time);
// myOutputWriterText.writeFile(current_time);

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
  
  // myOutputWriterParaview.writeFile(current_time);
  // myOutputWriterText.writeFile(current_time);
Iterationszahl=Iterationszahl+1;
}
std::cout<< "Noetige Iterationen " << Iterationszahl <<std::endl;
time1=clock()-tstart;
time1=time1/CLOCKS_PER_SEC;

std::cout<< "Laufzeit in s " << time1  <<std::endl;


}


