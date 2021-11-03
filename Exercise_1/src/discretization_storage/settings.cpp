#include "settings.h"
#include <fstream>
#include <iomanip>


void Settings::loadFromFile(std::string filename)
{
  float mA;
  // open file
  std::ifstream file(filename.c_str(), std::ios::in);

  // check if file is open
  if (!file.is_open())
  {
    std::cout << "Could not open parameter file \"" << filename << "\"." << std::endl;
    return;
  }
  // loop over lines of file
  while (!file.eof())
  {

    // read line

    std::string line;

    getline(file, line);

    std::string parameterName = "";
    std::string parameterValue = "";
    

    // erase the spaces, tabs and comments
    if (line.find("#")!=std::string::npos)
    {
      line.erase(line.find_first_of("#"));  
    }        
    while (line[0] == ' ' || line[0] == '\t')
    {
      line.erase(0,1);
    }

    int eqpos = line.find("=");
    if (eqpos<0)
    {
      continue;
    }
    
    parameterName =  line.substr(0,eqpos);
    parameterValue = line.substr(eqpos+1); 
    // print line if the parameter name is not empty
    if (parameterName != "")
    {
      int space=parameterName.find_first_of(" \t");
      if (space>0)
      {
      parameterName.erase(space);
      }
      while (parameterValue[0] == ' ' || parameterValue[0] == '\t')
    {
      parameterValue.erase(0,1);
    }
      space=parameterValue.find_first_of(" \t");
      if (space>0)
      {
      parameterValue.erase(space);
      }       
      //float pV=convertstringtonumber(parameterValue);
      //std::cout <<parameterName<< " with value " <<"<"<< pV << ">" <<std::endl;

      //std::cout <<parameterName<< " with value " <<"<"<< parameterValue << ">" <<std::endl;
       
      
      if (parameterName == "endTime")
      {
        endTime = atof(parameterValue.c_str());
      } 
      else if (parameterName =="physicalSizeX")
      {
        physicalSize[0]=atof(parameterValue.c_str());
      } 
      else if (parameterName =="physicalSizeY")
      {
        physicalSize[1] =atof(parameterValue.c_str());
      } 
      else if (parameterName =="re")
      {
        re =atof(parameterValue.c_str());
      } 
      else if (parameterName =="gX")
      {
        g[0] =atof(parameterValue.c_str());
      } 
      else if (parameterName =="gY")
      {
        g[1] =atof(parameterValue.c_str());
      } 
      else if (parameterName =="dirichletBottomX")
      {
        dirichletBcBottom[0] =atof(parameterValue.c_str());
      } 
      else if (parameterName =="dirichletBottomY")
      {
        dirichletBcBottom[1] =atof(parameterValue.c_str());
      } 
      else if (parameterName =="dirichletTopX")
      {
        dirichletBcTop[0] =atof(parameterValue.c_str());
      } 
      else if (parameterName =="dirichletTopY")
      {
        dirichletBcTop[1] =atof(parameterValue.c_str());
      } 
      else if (parameterName =="dirichletLeftX")
      {
        dirichletBcLeft[0] =atof(parameterValue.c_str());
      } 
      else if (parameterName =="dirichletLeftY")
      {
        dirichletBcLeft[1] =atof(parameterValue.c_str());
      } 
      else if (parameterName =="dirichletRightX")
      {
        dirichletBcRight[0] =atof(parameterValue.c_str());
      } 
      else if (parameterName =="dirichletRightY")
      {
        dirichletBcRight[1] =atof(parameterValue.c_str());
      } 
      else if (parameterName =="nCellsX")
      {
        nCells[0] =atoi(parameterValue.c_str());
      } 
      else if (parameterName =="nCellsY")
      {
        nCells[1] =atoi(parameterValue.c_str());
      }
      else if (parameterName =="useDonorCell")
      {
          if (parameterValue == "true")
          {
            useDonorCell=true;        
          } 
          else if (parameterValue == "false")
          {
            useDonorCell=false; 
          }
       }
       else if (parameterName =="alpha")
      {
        alpha =atof(parameterValue.c_str());
      } 
      else if (parameterName == "tau")
      {
        tau =atof(parameterValue.c_str());
      } 
      else if (parameterName =="maximumDt")
      {
        maximumDt =atof(parameterValue.c_str());
      } 
      else if (parameterName == "pressureSolver")
      {
      pressureSolver = parameterValue;
      } 
      else if (parameterName =="omega")
      {
        omega =atof(parameterValue.c_str());
      } 
      else if (parameterName =="epsilon")
      {
        epsilon =atof(parameterValue.c_str());
      } 
      else if (parameterName =="maximumNumberOfIterations")
      {
        // mit atoi kann er keine 1e umwandeln
        //mA=atof(parameterValue.c_str());
        //maximumNumberOfIterations = (int)mA;
        maximumNumberOfIterations=atoi(parameterValue.c_str());
      }
      
    
    }
  }
}


void Settings::printSettings()
{
  std::cout << "Settings: " << std::endl
    << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl
    << "  endTime: " << endTime << " s, re: " << re << ", g: (" << g[0] << "," << g[1] << "), tau: " << tau << ", maximum dt: " << maximumDt << std::endl
    << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1]  << ")"
    << ", top: ("  << dirichletBcTop[0] << "," << dirichletBcTop[1]  << ")"
    << ", left: ("  << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"
    << ", right: ("  << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << std::endl
    << "  useDonorCell: " << std::boolalpha << useDonorCell << ", alpha: " << alpha << std::endl
    << "  pressureSolver: " << pressureSolver << ", omega: " << omega << ", epsilon: " << epsilon << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
}

