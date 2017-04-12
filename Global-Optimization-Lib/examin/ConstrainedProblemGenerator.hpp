/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2016 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      ConstrainedProblemGenerator.h                               //
//                                                                         //
//  Purpose:   Header file for ExaMin problem interface                    //
//                                                                         //
//                                                                         //
//  Author(s): Lebedev I. Sovrasov V.                                      //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/**
\file ConstrainedProblemGenerator.h

\authors ������� �. �������� �.
\date 2016
\copyright ���� ��. �.�. ������������

\brief ����� ��� �������� ����� � �������������

\details
*/

#ifndef CONSTRAINED_PROBLEM_GENERATOR_H
#define CONSTRAINED_PROBLEM_GENERATOR_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <iostream>

#include "problemPar.hpp"
#include "ConstrainedProblem.hpp"

enum GenerateMode { RHS = 0, DELTA = 1 };

enum GeneratorOptions
{
  IMPROVE_OBJECTIVE = 1,
  TOTAL_DELTA = 1 << 1,
  ZOOM = 1 << 2,
  SHIFT = 1 << 3
};

template <class FType>
class TConstrainedProblemGenerator : virtual public TProblemPar
{
protected:

  /// ������� �������
  FType* mPObjective;
  /// �����������
  std::vector<FType*> mPConstraints;
  /// ��������� �����������: ���� ����� ���� ���� �������
  std::vector<double> mConstraintsParams;
  /// ����������� ������ ����������� ��� ���������� ������� �������
  std::vector<double> mImprovementOfTheObjective;
  /// ����� �� ��������� �����
  std::vector<bool> mNeedTuneParam;
  /// ����������� ��������������� ��� ������� �����������
  double** mvZoomRatios;
  /// �������������� ����� �����������
  double*** mvShift;
  /// ����� �� �������������� ������
  bool mIsZoom;
  /// ����� �� �������� ������� �����������
  bool mIsShift;

  /// ���� ����� ������(���� �������) ��� ���� ������� ��� ��� ������ ������� ����
  bool mIsTotalDelta;

  /// �������� �� ������� ������� ����� ���������� ����������� �� �����������
  bool mIsImprovementOfTheObjective;

  /// ������� ����������� �������
  FType* currentFunction;
  /// ����� ����������� �������
  int currentFunctionNumber;


  /// ��������� ���� ������� currentFunction
  double OneFunctionCalculate(double* y);
  /// �������� �������� ���� ����������� � �����
  double MaxFunctionCalculate(double* y);

  double FunctionCalculate(double* x);

  double CalculateRHS(double delta, int m = 100, double Epsilon = 0.01, int maxM = 10000000);

  /// ������ ����������� ��������������� ��� ������� �����������
  virtual void SetZoom();
  /// ������ ����� � ����������� �������� ��� ������� �����������
  virtual void SetShift();

  /** �������������� ������� � ������� index,
  �������������� ���� ��������� �������� ��������� ������� � ����������� �� ���������� ������
  */
  virtual void InitFunc(FType* func, int index)
  {  }

public:

  TConstrainedProblemGenerator();



  /** ����� ���������� ���������� ����� ����������� �������� ������� �������
  \param[out] y �����, � ������� ����������� ����������� ��������
  \return ��� ������ (#OK ��� #UNDEFINED)
  */
  virtual int GetOptimumPoint(double* y) const;

  /// ���������� ����� ����������� �������� ��� ������� fNumber
  virtual int GetConstraintOptimumPoint(double* point, int fNumber);

  /// �������� ����������� � ������
  void AddConstraint(FType* function, double parameter, int mode = DELTA, double imp = 0);

  /// ������ ������� �������
  void SetObjective(FType* function);

  /** ������� ������
  \param[in] generateOptions - ����, ���������� ����� ���������� �����. ��������
  ����������� ��������� ��������:
   GeneratorOptions::IMPROVE_OBJECTIVE - ������������ �� ����� ������� �������
   GeneratorOptions::TOTAL_DELTA - �������� ��� ����������� ������ ��� �� ����������� (���� true �� �������� mode � AddConstraint �� ������������)
   GeneratorOptions::ZOOM - �������������� �� ����������� ��� ��������� ������������� ������
   GeneratorOptions::SHIFT - �������� �� ����������� ��� �� ���������� ������� ��� ������ ���������� �������
  \return ����� ������
  */
  TConstrainedProblem<FType> GenerateProblem(int generateOptions);

};

/** ����� ���������� ���������� ����� ����������� �������� ������� �������
\param[out] y �����, � ������� ����������� ����������� ��������
\return ��� ������ (#OK ��� #UNDEFINED)
*/
template <class FType>
int TConstrainedProblemGenerator<FType>::GetOptimumPoint(double* y) const
{
  mPObjective->GetOptimumPoint(y);

  for (int i = 0; i < mDim; i++)
  {
    y[i] = y[i] * (*mvZoomRatios)[mNumberOfConstraints] + (*mvShift)[mNumberOfConstraints][i];
  }

  return TProblemPar::ProblemOK;
}

/// ���������� ����� ����������� �������� ��� ������� fNumber
template <class FType>
int TConstrainedProblemGenerator<FType>::GetConstraintOptimumPoint(double* point, int fNumber)
{
  mPConstraints[fNumber]->GetOptimumPoint(point);
  for (int i = 0; i < mDim; i++)
  {
    point[i] = point[i] * (*mvZoomRatios)[fNumber] + (*mvShift)[fNumber][i];
  }
  return TProblemPar::ProblemOK;
}

template <class FType>
TConstrainedProblemGenerator<FType>::TConstrainedProblemGenerator()
{
  mIsTotalDelta = true;
  mIsZoom = false;
  mIsShift = false;

  mvZoomRatios = 0;
  mvShift = 0;

  currentFunction = 0;
  currentFunctionNumber = 0;

  mPObjective = 0;
}

template <class FType>
void TConstrainedProblemGenerator<FType>::SetObjective(FType* function)
{
  mPObjective = function;
  mNumberOfCriterions = 1;
  mDim = function->GetDimension();
}

template <class FType>
double TConstrainedProblemGenerator<FType>::OneFunctionCalculate(double* y)
{
  double x[MAX_DIM];
  for (int i = 0; i < mDim; i++)
  {
    x[i] = (y[i] - (*mvShift)[currentFunctionNumber][i]) / (*mvZoomRatios)[currentFunctionNumber];
  }
  return currentFunction->Calculate(x);
}


/// �������� �������� ���� ����������� � �����
template <class FType>
double TConstrainedProblemGenerator<FType>::MaxFunctionCalculate(double* y)
{
  double x[MAX_DIM];
  if (GetRealNumberOfConstraints() > 0)
  {
    for (int j = 0; j < mDim; j++)
    {
      x[j] = (y[j] - (*mvShift)[0][j]) / (*mvZoomRatios)[0];
    }
    double res = mPConstraints[0]->Calculate(x);
    for (int i = 1; i < GetRealNumberOfConstraints(); i++)
    {
      for (int j = 0; j < mDim; j++)
      {
        x[j] = (y[j] - (*mvShift)[i][j]) / (*mvZoomRatios)[i];
      }

      double f = mPConstraints[i]->Calculate(x);
      if (res < f)
        res = f;
    }
    return res;
  }
  return 0;
}


template <class FType>
double TConstrainedProblemGenerator<FType>::FunctionCalculate(double* x)
{
  if (mIsTotalDelta)
    return MaxFunctionCalculate(x);
  else
    return OneFunctionCalculate(x);
}

template <class FType>
double TConstrainedProblemGenerator<FType>::CalculateRHS(double delta, int m, double Epsilon, int maxM)
{
  double rhs = 0;
  double* objectiveMin  = new double [mDim];
  unsigned dimension = mDim;
  double hmin = 0;

  for (int j = 0; j < mDim; j++)
    objectiveMin[j] = 0;
  GetOptimumPoint(objectiveMin);

  hmin = mPObjective->GetOptimumValue();
  //mPObjective->GetOptimumValue(hmin);

  double hmax = hmin;
  double d = 0;
  //����������� �������, � ����� - ���������
  int* size = new int[dimension];//���-�� ����� ����� �� �����������
  double* step = new double[dimension];//��� �� ������ �����������
  int sumn = 1;//����� ���������

  double* a = new double[dimension];
  double* b = new double[dimension];
  mPObjective->GetBounds(a, b);
  double multiplyingLength = 1;
  for (unsigned i = 0; i < dimension; i++)
  {
    d = (b[i] - a[i]);
    size[i] = (int)ceil(d / Epsilon) + 1;
    step[i] = d / (size[i] - 1);
    multiplyingLength = multiplyingLength * d;
    sumn *= (size[i]);
  }

  if ((sumn > maxM) || (sumn <= 0))
  {
    multiplyingLength = multiplyingLength / maxM;
    Epsilon = pow(multiplyingLength, 1.0 / (double) dimension);
    sumn = 1;
    multiplyingLength = 1;

    for (unsigned i = 0; i < dimension; i++)
    {
      d = (b[i] - a[i]);
      size[i] = (int)ceil(d / Epsilon) + 1;
      step[i] = d / (size[i] - 1);
      sumn *= (size[i]);
    }

  }

  double* f = new double[sumn];//�������� �������
  double* yArray = new double[dimension*omp_get_max_threads()];

  #pragma omp parallel for num_threads(omp_get_max_threads())
  for (int i = 0; i < sumn; i++)
  {
    double w;
    int z = i;
    double* y = yArray + omp_get_thread_num()*dimension;
    //��������� ���������� ����� ���������
    for (unsigned j = 0; j < dimension; j++)
    {
      w = z % size[j];//���������� ����� ���� �� i-�� �����������
      y[j] = a[j] + w * step[j];//����� ������� + ����� ���� �� i-�� ����������� * ��� �� i-�� �����������
      z = z / size[j];//��� ���������� ������ ���� �� ��������� �����������
    }
    //�������� ���������
    f[i] = FunctionCalculate(y);
    hmax = std::max(f[i], hmax);
    hmin = std::min(f[i], hmin);
  }

  double* h1 = new double[m];
  double* h2 = new double[m];
  int* p = new int[m];
  int* s = new int[m];

  double deltah = (hmax - hmin) / m;

  for (int i = 0; i < m; i++)
  {
    h1[i] = hmin + i * deltah;
    h2[i] = hmin + (i + 1) * deltah;
    p[i] = 0;
    s[i] = 0;
  }

  for (int i = 0; i < sumn; i++)
    for (int j = 0; j < m; j++)
      if ((f[i] >= h1[j]) && (f[i] <= h2[j]))
      {
        p[j] ++;
        break;
      }

      s[0] = p[0];
      for (int i = 1; i < m; i++)
      {
        s[i] = s[i - 1] + p[i];
      }

      double smax = s[m - 1];
      double g = delta * smax;
      for (int i = 0; i < m; i++)
      {
        if (s[i] >= g)
        {
          rhs = h2[i];
          break;
        }
      }

      double dm = delta;
      if (dm == 0)
        dm += 0.1;
      dm = dm * (hmax - hmin);

      double criticalValue = FunctionCalculate(objectiveMin);

      if (rhs < criticalValue)
      {
        std::cout << "Q was changed from " << rhs << " to " << criticalValue + dm << "\n";
        rhs = criticalValue + dm;
      }

      delete[] size;
      delete[] step;
      delete[] a;
      delete[] b;
      delete[] f;
      delete[] yArray;

      delete[] h1;
      delete[] h2;
      delete[] p;
      delete[] s;

      return rhs;
}

template <class FType>
void TConstrainedProblemGenerator<FType>::AddConstraint(FType* function, double parameter, int mode, double imp)
{
  mPConstraints.push_back(function);
  mConstraintsParams.push_back(parameter);
  if (mode == DELTA)
    mNeedTuneParam.push_back(true);
  else
    mNeedTuneParam.push_back(false);
  mImprovementOfTheObjective.push_back(imp);
  mNumberOfConstraints++;
}

template <class FType>
TConstrainedProblem<FType> TConstrainedProblemGenerator<FType>::
  GenerateProblem(int generateOptions)
{
  mIsShift = generateOptions & SHIFT;
  mIsZoom = generateOptions & ZOOM;
  /// �������, � ������ �����������, ����� �������
  FType** mPFunction = new FType* [mNumberOfConstraints + 1];
  for (int i = 0; i < mNumberOfConstraints; i++)
  {
    mPFunction[i] = mPConstraints[i];
  }
  mPFunction[mNumberOfConstraints] = mPObjective;

  double* objectiveMin = new double [mDim];

  /// ����� �����������, ��� �� RHS
  double** mQ = new double* [1];
  (*mQ) = new double [mNumberOfConstraints + 1];

  /// ����������� ��������������� ��� ������� �����������
  double* mZoomRatios = new double [mNumberOfConstraints + 1];
  /// �������������� ����� �����������
  double** mShift = new double* [mNumberOfConstraints + 1];

  for (int i = 0; i < GetRealNumberOfFunctions(); i++)
  {
    mZoomRatios[i] = 1.0;

    mShift[i] = new double [mDim];
    for (int j = 0; j < mDim; j++)
    {
      mShift[i][j] = 0.0;
    }
  }

  mvZoomRatios = &mZoomRatios;
  mvShift = &mShift;

  mIsImprovementOfTheObjective = generateOptions & IMPROVE_OBJECTIVE;
  mIsTotalDelta = generateOptions & TOTAL_DELTA;
  /// ����������� ���������
  double** improvementCoefficients = new double* [1];
  *improvementCoefficients = new double [mNumberOfConstraints + 1];
  for (int i = 0; i < mNumberOfConstraints; i++)
    (*improvementCoefficients)[i] = mImprovementOfTheObjective[i];
  (*improvementCoefficients)[mNumberOfConstraints] = 0;

  for (int i = 0; i < GetRealNumberOfFunctions(); i++)
  {
    InitFunc(mPFunction[i], i);
  }

  if (mIsZoom)
  {
    SetZoom();
  }

  if (mIsShift)
  {
    SetShift();
  }

  currentFunctionNumber = 0;

  if (mIsTotalDelta)
  {
    double q = CalculateRHS(mConstraintsParams[0]);
    for (int i = 0; i < GetRealNumberOfConstraints(); i++)
    {
      (*mQ)[i] = q;
    }
  }
  else
  {
    for (int i = 0; i < GetRealNumberOfConstraints(); i++)
    {
      if (mNeedTuneParam[i])
      {
        currentFunctionNumber = i;
        currentFunction = mPFunction[i];
        (*mQ)[i] = CalculateRHS(mConstraintsParams[i]);
      }
      else
      {
        (*mQ)[i] = mConstraintsParams[i];
      }

    }
  }

  return TConstrainedProblem<FType>(mPFunction, mQ, mZoomRatios, mShift,
    mIsImprovementOfTheObjective, improvementCoefficients, mNumberOfConstraints, mNumberOfCriterions, mDim);
}

/// ������ ����������� ��������������� ��� ������� �����������
template <class FType>
void TConstrainedProblemGenerator<FType>::SetZoom()
{

  double* objectiveMin = new double [mDim];
  double* lower = new double [mDim];
  double* upper = new double [mDim];
  double* constraintMin = new double [mDim];

  for (int j = 0; j < mDim; j++)
  {
    objectiveMin[j] = 0;
    lower[j] = 0;
    upper[j] = 0;
  }
  /// ���������� ���������� ����������� �������� ������� �������
  GetOptimumPoint(objectiveMin);
  /// ���������� ������� �������
  mPObjective->GetBounds(lower, upper);

  double maxDistanceToBoundary = 0;

  for (int k = 0; k < mDim; k++)
  {
    if (maxDistanceToBoundary < (objectiveMin[k] - lower[k]))
      maxDistanceToBoundary = (objectiveMin[k] - lower[k]);
    if (maxDistanceToBoundary < (upper[k] - objectiveMin[k]))
      maxDistanceToBoundary = (upper[k] - objectiveMin[k]);
  }

  if (fabs(maxDistanceToBoundary) < AccuracyDouble )
  {
    mIsZoom = false;
    for (int j = 0; j < GetRealNumberOfFunctions(); j++)
    {
      (*mvZoomRatios)[j] = 1;
    }

    delete [] objectiveMin;
    delete [] lower;
    delete [] upper;
    delete [] constraintMin;
    return;
  }

  for (int i = 0; i < GetRealNumberOfConstraints(); i++)
  {
    double minDistanceToBoundary = upper[0] - lower[0];

    for (int j = 0; j < mDim; j++)
    {
      constraintMin[j] = 0;
    }

    GetConstraintOptimumPoint(constraintMin, i);
    for (int k = 0; k < mDim; k++)
    {
      if (minDistanceToBoundary > (constraintMin[k] - lower[k]))
        minDistanceToBoundary = (constraintMin[k] - lower[k]);
      if (minDistanceToBoundary > (upper[k] - constraintMin[k]))
        minDistanceToBoundary = (upper[k] - constraintMin[k]);
    }

    if (fabs(minDistanceToBoundary) < AccuracyDouble )
    {
      mIsZoom = false;
      for (int j = 0; j < GetRealNumberOfFunctions(); j++)
      {
        (*mvZoomRatios)[j] = 1;
      }
      delete [] objectiveMin;
      delete [] lower;
      delete [] upper;
      delete [] constraintMin;
      return;
    }
    else
    {
      (*mvZoomRatios)[i] = maxDistanceToBoundary / minDistanceToBoundary;
    }
  }

  delete [] objectiveMin;
  delete [] lower;
  delete [] upper;
  delete [] constraintMin;
}

/// ������ ����� � ����������� �������� ��� ������� �����������
template <class FType>
void TConstrainedProblemGenerator<FType>::SetShift()
{
  double* objectiveMin = new double [mDim];
  double* lower = new double [mDim];
  double* upper = new double [mDim];
  double* constraintMin = new double [mDim];

  for (int j = 0; j < mDim; j++)
  {
    objectiveMin[j] = 0;
    lower[j] = 0;
    upper[j] = 0;
  }
  /// ���������� ���������� ����������� �������� ������� �������
  GetOptimumPoint(objectiveMin);
  /// ���������� ������� �������
  mPObjective->GetBounds(lower, upper);

  for (int i = 0; i < GetRealNumberOfConstraints(); i++)
  {

    for (int j = 0; j < mDim; j++)
    {
      constraintMin[j] = 0;
    }

    GetConstraintOptimumPoint(constraintMin, i);

    for (int k = 0; k < mDim; k++)
    {
      (*mvShift)[i][k] = objectiveMin[k] - constraintMin[k];
    }

    GetConstraintOptimumPoint(constraintMin, i);
  }

  delete [] objectiveMin;
  delete [] lower;
  delete [] upper;
  delete [] constraintMin;
}

#endif // CONSTRAINED_PROBLEM_GENERATOR_H
