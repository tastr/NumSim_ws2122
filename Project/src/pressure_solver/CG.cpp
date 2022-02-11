#include "CG.h"
#include <cassert>

CG::CG(Discretization &discretization_)
    : PressureSolver(discretization_),
      Matrix({discretization_.getNumberOfFluidCerlls(), discretization_.getNumberOfFluidCerlls()}),
      pvektor({discretization_.getNumberOfFluidCerlls(), 1}),
      rhsVektor({discretization_.getNumberOfFluidCerlls(), 1})
{
    setMatrix();
    // Array2D k = matMulVec(Matrix, pvektor);
    // matrixPrint(k);
}

CG::~CG()
{
}

void CG::calculateP()
{
    int i;
    int j;

    double epsilonquad = discretization_.getepsilon() * discretization_.getepsilon();
    setRHSVektor();
    int iter = 0;
    double sum;
    double sum2;
    Array2D r0(pvektor.size());
    Array2D r1(pvektor.size());
    Array2D d0(pvektor.size());
    Array2D d1(pvektor.size());
    Array2D z(pvektor.size());
    Array2D pnew(pvektor.size());
    double a = 0;
    double b = 0;
    double resterm;

    
    do
    {
        for (int j = 0; j < Matrix.size()[0]; j++)
        {
            pnew(j, 0) = rhsVektor(j, 0);
        }

        for (int j = 0; j < Matrix.size()[0]; j++)
        {
            for (int i = 0; i < Matrix.size()[0]; i++)
            {
                if (j != i)
                {
                    pnew(j, 0) = pnew(j, 0) - Matrix(i, j) * pvektor(i, 0);
                }
            }
            pnew(j, 0) = pnew(j, 0) / Matrix(j, j);//from D^-1


        }

        for (int j = 0; j < Matrix.size()[0]; j++)
        {
            pvektor(j, 0) = pnew(j, 0);
        }

        for (int j = 0; j < pvektor.size()[0]; j++)
        {
            sum = 0;
            for (int i = 0; i < Matrix.size()[0]; i++)
            {
                sum += pvektor(i, 0) * Matrix(i, j);
            }
            r1(j, 0) = rhsVektor(j, 0) - sum;
        }

        iter += 1;
      // printf("Norm vec %10.8f \n", normVec(r1));

       }while (normVec(r1) > epsilonquad && iter < discretization_.getMaxIteration());
   //  printf("Norm vec %f tol %f Iter %d  \n", normVec(r1),epsilonquad,iter);
  std::cout << "Residuum " << normVec(r1)<<" tol " << epsilonquad<< " Safe " << iter << std::endl;
  for (int n = 0; n < Matrix.size()[0]; n++)
        {
            i = discretization_.getFluidCellsIndices(n)[0];
            j = discretization_.getFluidCellsIndices(n)[1];

            discretization_.setP(i, j, pvektor(n, 0));
        }
        discretization_.updatedPressureBC();
        discretization_.setObstaclePressure();
   
    
  // resterm = residuum() / Matrix.size()[0];
  //std::cout << "Residuum " << residuum() << " Safe " << iter << std::endl;

  // }while (resterm > epsilonquad && iter < 20000);
     
    //std::cout << "Residuum " << residuum() << " Safe " << iter << std::endl;


/* for (int j = 0; j < pvektor.size()[0]; j++)
{
    sum = 0;
    for (int i = 0; i < Matrix.size()[0]; i++)
    {
        sum += pvektor(i, 0) * Matrix(i, j);
    }
    r0(j,0)=rhsVektor(j,0)- sum;
    d0(j,0)=r0(j,0);
}
do
{

  for (int j = 0; j < pvektor.size()[0]; j++)
{
    sum = 0;
    for (int i = 0; i < Matrix.size()[0]; i++)
    {
        sum += d0(i, 0) * Matrix(i, j);
    }
    z(j,0)=sum;
}
sum=0;
sum2=0;

 for (int j = 0; j < pvektor.size()[0]; j++)
{
 sum+=r0(j,0)*r0(j,0);
 sum2+=r0(j,0)*d0(j,0);
}
a=sum/sum2;

 for (int j = 0; j < pvektor.size()[0]; j++)
{
    pvektor(j,0)+=a*d0(j,0);
    r1(j,0)=r0(j,0)-a*z(j,0);
}

sum=0;
sum2=0;

 for (int j = 0; j < pvektor.size()[0]; j++)
{
 sum+=r1(j,0)*r1(j,0);
 sum2+=r0(j,0)*r0(j,0);
}
b=sum/sum2;
 for (int j = 0; j < pvektor.size()[0]; j++)
{
 d0(j,0)=b*d0(j,0)+r1(j,0);
 r0(j,0)=r1(j,0);
}


iter+=1;
 printf("%10.5f \n",normVec(r1));
} while (normVec(r1)>epsilonquad && iter< discretization_.getMaxIteration());

*/
}

void CG::setMatrix()
{
    int i, j, k;
    double dy = discretization_.dy();
    double dx = discretization_.dx();
    double dy2 = discretization_.dy() * discretization_.dy();
    double dx2 = discretization_.dx() * discretization_.dx();

    for (int n = 0; n < Matrix.size()[0]; n++)
    {
        i = discretization_.getFluidCellsIndices(n)[0];
        j = discretization_.getFluidCellsIndices(n)[1];
        printf("i=%d j=%d\n", i, j);
        Matrix(n, n) = -2 * (1 / dx2 + 1 / dy2);

        // Check if its a border Point
        if (isBorder(i, j))
        {
            if (isEdge(i, j))
            {
                setEdge(i, j, n);
            }
            else
            {
                setBorder(i, j, n);
            }
        }
        else
        {
            setInnerValue(i, j, n);
        }
    }
    printf("%d  %d \n", Matrix.size()[0], Matrix.size()[1]);
    printf("%d \n", discretization_.nCells()[1]);
}

bool CG::isEdge(int i, int j)
{
    bool edge = false;
    if (i == discretization_.nCells()[0] || i == 1)
    {
        if (j == discretization_.nCells()[1] || j == 1)
        {
            edge = true;
        }
    }
    return edge;
}

bool CG::isBorder(int i, int j)
{
    bool border = false;
    if (i == discretization_.nCells()[0] || i == 1 || j == discretization_.nCells()[1] || j == 1)
    {
        border = true;
    }
    return border;
}

void CG::setEdge(int i, int j, int n)
{
    if (i == 1 && j == 1)
    {
        setLowerLeftEdge(i, j, n);
    }
    else if (i == discretization_.nCells()[1] && j == 1)
    {
        setLowerRightEdge(i, j, n);
    }
    else if (i == 1 && j == discretization_.nCells()[1])
    {
        setUpperLeftEdge(i, j, n);
    }
    else if (i == discretization_.nCells()[0] && j == discretization_.nCells()[1])
    {
        setUpperRightEdge(i, j, n);
    }
}

void CG::setBorder(int i, int j, int n)
{
    double dy2 = discretization_.dy() * discretization_.dy();
    double dx2 = discretization_.dx() * discretization_.dx();

    if (i == 1)
    {
        if (!discretization_.isOutFlowLeft())
        {
            Matrix(n, n) += 1 / dx2;
        }
        setValueTop(i, j, n);
        setValueBottom(i, j, n);
        setValueRight(i, j, n);
    }
    else if (i == discretization_.nCells()[0])
    {
        if (!discretization_.isOutFlowRight())
        {
            Matrix(n, n) += 1 / dx2;
        }
        setValueTop(i, j, n);
        setValueBottom(i, j, n);
        setValueLeft(i, j, n);
    }

    if (j == 1)
    {
        if (!discretization_.isOutFlowBottom())
        {
            Matrix(n, n) += 1 / dy2;
        }
        setValueTop(i, j, n);
        setValueLeft(i, j, n);
        setValueRight(i, j, n);
    }
    else if (j == discretization_.nCells()[1])
    {
        if (!discretization_.isOutFlowTop())
        {
            Matrix(n, n) += 1 / dy2;
        }
        setValueRight(i, j, n);
        setValueBottom(i, j, n);
        setValueLeft(i, j, n);
    }
}

void CG::setInnerValue(int i, int j, int n)
{
    setValueRight(i, j, n);
    setValueBottom(i, j, n);
    setValueLeft(i, j, n);
    setValueTop(i, j, n);
}

void CG::setValueTop(int i, int j, int n)
{
    double dy2 = discretization_.dy() * discretization_.dy();
    if (discretization_.getTyp(i, j + 1) == 0)
    {
        Matrix(discretization_.getFluidcellsIndex(i, j + 1), n) = 1/dy2;
    }
    else
    {
        Matrix(n, n) += 1 / dy2;
    }
}
void CG::setValueLeft(int i, int j, int n)
{
    double dx2 = discretization_.dx() * discretization_.dx();
    if (discretization_.getTyp(i - 1, j) == 0)
    {
        Matrix(discretization_.getFluidcellsIndex(i - 1, j), n) = 1/dx2;
    }
    else
    {
        Matrix(n, n) += 1 / dx2;
    }
}
void CG::setValueBottom(int i, int j, int n)
{
    double dy2 = discretization_.dy() * discretization_.dy();
    if (discretization_.getTyp(i, j - 1) == 0)
    {
        Matrix(discretization_.getFluidcellsIndex(i, j - 1), n) = 1/dy2;
    }
    else
    {
        Matrix(n, n) += 1 / dy2;
    }
}
void CG::setValueRight(int i, int j, int n)
{
    double dx2 = discretization_.dx() * discretization_.dx();
    if (discretization_.getTyp(i + 1, j) == 0)
    {
        Matrix(discretization_.getFluidcellsIndex(i + 1, j), n) = 1/dx2;
    }
    else
    {
        Matrix(n, n) += 1 / dx2;
    }
}

void CG::setUpperLeftEdge(int i, int j, int n)
{
    double dy2 = discretization_.dy() * discretization_.dy();
    double dx2 = discretization_.dx() * discretization_.dx();

    if (!discretization_.isOutFlowLeft() && !discretization_.isOutFlowTop())
    {
        Matrix(n, n) += 1 / dx2 + 1 / dy2;
    }
    else if (!discretization_.isOutFlowLeft())
    {
        Matrix(n, n) += 1 / dy2;
    }
    else if (!discretization_.isOutFlowTop())
    {
        Matrix(n, n) += 1 / dx2;
    }
    setValueRight(i, j, n);
    setValueBottom(i, j, n);
}

void CG::setUpperRightEdge(int i, int j, int n)
{
    double dy2 = discretization_.dy() * discretization_.dy();
    double dx2 = discretization_.dx() * discretization_.dx();

    if (!discretization_.isOutFlowRight() && !discretization_.isOutFlowTop())
    {
        Matrix(n, n) += 1 / dx2 + 1 / dy2;
    }
    else if (!discretization_.isOutFlowRight())
    {
        Matrix(n, n) += 1 / dy2;
    }
    else if (!discretization_.isOutFlowTop())
    {
        Matrix(n, n) += 1 / dx2;
    }
    setValueLeft(i, j, n);
    setValueBottom(i, j, n);
}

void CG::setLowerLeftEdge(int i, int j, int n)
{
    double dy2 = discretization_.dy() * discretization_.dy();
    double dx2 = discretization_.dx() * discretization_.dx();
    if (!discretization_.isOutFlowLeft() && !discretization_.isOutFlowBottom())
    {
        Matrix(n, n) += 1 / dx2 + 1 / dy2;
    }
    else if (!discretization_.isOutFlowLeft())
    {
        Matrix(n, n) += 1 / dy2;
    }
    else if (!discretization_.isOutFlowBottom())
    {
        Matrix(n, n) += 1 / dx2;
    }
    setValueRight(i, j, n);
    setValueTop(i, j, n);
}

void CG::setLowerRightEdge(int i, int j, int n)
{
    double dy2 = discretization_.dy() * discretization_.dy();
    double dx2 = discretization_.dx() * discretization_.dx();

    if (!discretization_.isOutFlowRight() && !discretization_.isOutFlowBottom())
    {
        Matrix(n, n) += 1 / dx2 + 1 / dy2;
    }
    else if (!discretization_.isOutFlowRight())
    {
        Matrix(n, n) += 1 / dy2;
    }
    else if (!discretization_.isOutFlowBottom())
    {
        Matrix(n, n) += 1 / dx2;
    }
    setValueLeft(i, j, n);
    setValueTop(i, j, n);
}

void CG::setRHSVektor()
{
    int i, j;
    for (int n = 0; n < Matrix.size()[0]; n++)
    {
        i = discretization_.getFluidCellsIndices(n)[0];
        j = discretization_.getFluidCellsIndices(n)[1];
        rhsVektor(n, 0) = discretization_.rhs(i, j);
    }
}

Array2D CG::matMul(Array2D A, Array2D B)
{
    std::array<int, 2> sizeA = A.size();
    std::array<int, 2> sizeB = B.size();
    double sum;
    if (sizeA[0] != sizeB[1])
    {
        printf("Falsche Matrixdimension");
        return Array2D({1, 1});
    }
    Array2D result({sizeB[0], sizeA[1]});

    for (int j = 0; j < sizeB[1]; j++)
    {
        sum = 0;
        for (int k = 0; k < sizeB[0]; k++)
        {
            for (int i = 0; i < sizeA[0]; i++)
            {
                sum += B(k, j) * A(i, j);
            }
            result(k, j) = sum;
        }
    }
    return result;
}

Array2D CG::matMulVec(Array2D A, Array2D B)
{
    std::array<int, 2> sizeA = A.size();
    std::array<int, 2> sizeB = B.size();
    double sum;
    if (sizeA[0] != sizeB[0])
    {
        printf("Falsche Matrixdimension");
        return Array2D({1, 1});
    }
    Array2D result({sizeA[1], 1});

    for (int j = 0; j < sizeB[0]; j++)
    {
        sum = 0;
        for (int i = 0; i < sizeA[0]; i++)
        {
            sum += B(j, 0) * A(i, j);
        }
        result(j, 0) = sum;
    }
    return result;
}

double CG::matMulscal(Array2D A, Array2D B)
{
    std::array<int, 2> sizeA = A.size();
    std::array<int, 2> sizeB = B.size();
    double sum;
    if (sizeA[0] != sizeB[0])
    {
        printf("Falsche Matrixdimension");
        return 0;
    }
    double result = 0;

    for (int j = 0; j < sizeB[0]; j++)
    {
        result += A(j, 0) * B(j, 0);
    }
    return result;
}

Array2D CG::vecAdd(Array2D A, Array2D B)
{
    Array2D result(B.size());
    for (int j = 0; j < B.size()[0]; j++)
    {
        result(j, 0) = A(j, 0) + B(j, 0);
    }
    return result;
}

Array2D CG::vecSub(Array2D A, Array2D B)
{
    Array2D result(B.size());
    for (int j = 0; j < B.size()[0]; j++)
    {
        result(j, 0) = A(j, 0) - B(j, 0);
    }
    return result;
}

Array2D CG::vecMultScal(double a, Array2D A)
{
    Array2D result(A.size());
    for (int j = 0; j < A.size()[0]; j++)
    {
        result(j, 0) = a * A(j, 0);
    }
    return result;
}

void CG::matrixPrint(Array2D Matrix)
{
    for (int j = 0; j < Matrix.size()[1]; j++)
    {
        for (int i = 0; i < Matrix.size()[0]; i++)
        {
            printf("%3.0f ", Matrix(i, j));
        }
        printf("\n");
    }
}

double CG::normVec(Array2D A)
{
    double result = 0;
    for (int j = 0; j < A.size()[0]; j++)
    {
        result += A(j, 0) * A(j, 0);
    }
    return result / A.size()[0];
}
