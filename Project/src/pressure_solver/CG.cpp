#include "CG.h"
#include <cassert>

CG::CG(Discretization &discretization_)
    : PressureSolver(discretization_),
      Matrix({discretization_.getNumberOfFluidCerlls(), discretization_.getNumberOfFluidCerlls()}),
      pvektor({discretization_.getNumberOfFluidCerlls(),1}),
      rhsVektor({discretization_.getNumberOfFluidCerlls(),1})      
{
    setMatrix();
    Array2D k= matMulVec(Matrix, pvektor);
    matrixPrint(k);
}

CG::~CG()
{
}

void CG::calculateP()
{
setRHSVektor();
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
        Matrix(n, n) = -2 * (1 / (dx2) + 1 / dy2);

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
        Matrix(discretization_.getFluidcellsIndex(i, j + 1), n) = 1;
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
        Matrix(discretization_.getFluidcellsIndex(i - 1, j), n) = 1;
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
        Matrix(discretization_.getFluidcellsIndex(i, j - 1), n) = 1;
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
        Matrix(discretization_.getFluidcellsIndex(i + 1, j), n) = 1;
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
{int i,j;
    for (int n = 0; n < Matrix.size()[0]; n++)
    {
        i = discretization_.getFluidCellsIndices(n)[0];
        j = discretization_.getFluidCellsIndices(n)[1];
        rhsVektor(1,n)=discretization_.rhs(i,j);
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
    Array2D result({sizeB[0],sizeA[1]});

    for (int j = 0; j < sizeB[1]; j++)
    {
        sum=0;
        for (int k = 0; k < sizeB[0]; k++)
        {
            for (int i = 0; i < sizeA[0]; i++)
            {
             sum+=B(k,j)*A(i,j); 
            }
            result(k,j)=sum;
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
    Array2D result({sizeA[1],1});

    for (int j = 0; j < sizeB[0]; j++)
    {
        sum=0;
             for (int i = 0; i < sizeA[0]; i++)
            {
             sum+=B(j,0)*A(i,j); 
            }
            result(j,0)=sum;
        
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
    double result=0;

    for (int j = 0; j < sizeB[0]; j++)
    {
        result+=A(j,0)*B(j,0);
    }
    return result;
}

Array2D CG::vecAdd(Array2D A, Array2D B)
{ Array2D result(B.size());   
    for (int j = 0; j < B.size()[0]; j++)
    {
        result(j,0)=A(j,0)+B(j,0);
    }
    return result;
}

Array2D CG::vecMultScal(double a, Array2D A)
{ Array2D result(A.size());   
    for (int j = 0; j < A.size()[0]; j++)
    {
        result(j,0)=a*A(j,0);
    }
    return result;

}


void CG::matrixPrint(Array2D Matrix)
{    
    for (int j = 0; j < Matrix.size()[1]; j++){
    for (int i = 0; i < Matrix.size()[0]; i++)
    {
        printf("%3.0f ",Matrix(i,j));
    }
   printf("\n");
   }
}