/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2016 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      ConstrainedProblem.h                                        //
//                                                                         //
//  Purpose:   Header file for ExaMin problem interface                    //
//                                                                         //
//                                                                         //
//  Author(s): Lebedev I. Sovrasov V.                                      //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/**
\file ConstrainedProblem.h

\authors ������� �. �������� �.
\date 2016
\copyright ���� ��. �.�. ������������

\brief ����� ��� ���������� ����� � �������������

\details
*/

#ifndef CONSTRAINED_PROBLEM_H
#define CONSTRAINED_PROBLEM_H

#include <vector>
#include <cmath>
#include <algorithm>

#include "problemPar.hpp"

template <class FType>
class TConstrainedProblemBase : virtual public TProblemPar
{
protected:
  /// �������, � ������ �����������, ����� �������
  FType** mPFunction;
  /// ����� �����������, ��� �� RHS
  double** mQ;

  /// ����������� ��������������� ��� ������� �����������
  double* mZoomRatios;
  /// �������������� ����� �����������
  double** mShift;
  /// �������� �� ������� ������� ����� ���������� ����������� �� �����������
  bool mIsImprovementOfTheObjective;
  /// ����������� ���������
  double** mImprovementCoefficients;

public:

  TConstrainedProblemBase();

  TConstrainedProblemBase(const TConstrainedProblemBase& problem);

  TConstrainedProblemBase(FType** PFunction, double** q, double* ZoomRatios, double** Shift,
    bool IsImprovementOfTheObjective, double** ImprovementCoefficients,
    int CountConstrained, int NumberOfCriterions, int Dim);

  /// ���������� ����� �����������
  virtual double GetFunctionRHS(int fNumber) const;

  /// ����������� ������
  virtual int GetDimension() const;

  /// ���������� ������� ������� ������
  virtual void GetBounds(double *lb, double* ub);

  /// ���������� ����� �����������
  virtual int GetConstraintsNumber() const;

  /** ��������� ������� � ������� fNumber � ����� y
  � ������ �����������, ����� �������
  */
  virtual double CalculateFunction(const double* y, int fNumber);

};

template <class FType>
TConstrainedProblemBase<FType>::TConstrainedProblemBase()
{
};

template <class FType>
TConstrainedProblemBase<FType>::TConstrainedProblemBase(const TConstrainedProblemBase& problem)
{
  mNumberOfConstraints = problem.mNumberOfConstraints;
  mNumberOfCriterions = problem.mNumberOfCriterions;
  mDim = problem.mDim;
  mIsImprovementOfTheObjective = problem.mIsImprovementOfTheObjective;

  mPFunction = new FType* [mNumberOfConstraints + 1];
  for (int i = 0; i < mNumberOfConstraints + 1; i++)
    mPFunction[i] = problem.mPFunction[i];
  mQ = new double* [1];
  (*mQ) = new double [mNumberOfConstraints + 1];
  for (int i = 0; i < mNumberOfConstraints + 1; i++)
    (*mQ)[i] = (*problem.mQ)[i];

  mImprovementCoefficients = new double* [1];
  *mImprovementCoefficients = new double [mNumberOfConstraints + 1];
  for (int i = 0; i < mNumberOfConstraints + 1; i++)
    (*mImprovementCoefficients)[i] = (*problem.mImprovementCoefficients)[i];

  /// ����������� ��������������� ��� ������� �����������
  mZoomRatios = new double [mNumberOfConstraints + 1];
  /// �������������� ����� �����������
  mShift = new double* [mNumberOfConstraints + 1];

  for (int i = 0; i < GetRealNumberOfFunctions(); i++)
  {
    mZoomRatios[i] = problem.mZoomRatios[i];

    mShift[i] = new double [mDim];
    for (int j = 0; j < mDim; j++)
    {
      mShift[i][j] = problem.mShift[i][j];
    }
  }
}

template <class FType>
TConstrainedProblemBase<FType>::TConstrainedProblemBase(FType** PFunction, double** q, double* ZoomRatios, double** Shift,
  bool IsImprovementOfTheObjective, double** ImprovementCoefficients,
  int CountConstrained, int NumberOfCriterions, int Dim)
{
  mNumberOfConstraints = CountConstrained;
  mNumberOfCriterions = NumberOfCriterions;
  mDim = Dim;
  mIsImprovementOfTheObjective = IsImprovementOfTheObjective;

  mPFunction = new FType* [mNumberOfConstraints + 1];
  for (int i = 0; i < mNumberOfConstraints + 1; i++)
    mPFunction[i] = PFunction[i];
  mQ = new double* [1];
  (*mQ) = new double [mNumberOfConstraints + 1];
  for (int i = 0; i < mNumberOfConstraints + 1; i++)
    (*mQ)[i] = (*q)[i];

  mImprovementCoefficients = new double* [1];
  *mImprovementCoefficients = new double [mNumberOfConstraints + 1];
  for (int i = 0; i < mNumberOfConstraints + 1; i++)
    (*mImprovementCoefficients)[i] = (*ImprovementCoefficients)[i];

  /// ����������� ��������������� ��� ������� �����������
  mZoomRatios = new double [mNumberOfConstraints + 1];
  /// �������������� ����� �����������
  mShift = new double* [mNumberOfConstraints + 1];

  for (int i = 0; i < GetRealNumberOfFunctions(); i++)
  {
    mZoomRatios[i] = ZoomRatios[i];

    mShift[i] = new double [mDim];
    for (int j = 0; j < mDim; j++)
    {
      mShift[i][j] = Shift[i][j];
    }
  }
}

template <class FType>
double TConstrainedProblemBase<FType>::GetFunctionRHS(int fNumber) const
{
  return (*mQ)[fNumber];
}

template <class FType>
int TConstrainedProblemBase<FType>::GetDimension() const
{
  return mPFunction[mNumberOfConstraints]->GetDimension();
}

template <class FType>
void TConstrainedProblemBase<FType>::GetBounds(double *lb, double* ub)
{
  mPFunction[mNumberOfConstraints]->GetBounds(lb, ub);
}

template <class FType>
int TConstrainedProblemBase<FType>::GetConstraintsNumber() const
{
  return mNumberOfConstraints;
}

template <class FType>
double TConstrainedProblemBase<FType>::CalculateFunction(const double* y, int fNumber)
{
  double x[MAX_DIM];
  if (fNumber < GetConstraintsNumber())
  {
    for (int i = 0; i < mDim; i++)
    {
      x[i] = (y[i] - mShift[fNumber][i]) / mZoomRatios[fNumber];
    }
    return mPFunction[fNumber]->Calculate(x) - (*mQ)[fNumber];
  }
  else
  {
    if (mIsImprovementOfTheObjective)
    {
      double resultCoefficient = 0;

      for (int j = 0; j < GetRealNumberOfConstraints(); j++)
      {
        for (int i = 0; i < mDim; i++)
        {
          x[i] = (y[i] - mShift[j][i]) / mZoomRatios[j];
        }
        double f = mPFunction[j]->Calculate(x) - (*mQ)[j];
        double fVal = std::max(f, 0.0);
        resultCoefficient += (*mImprovementCoefficients)[j] * (fVal * fVal * fVal);
      }
      double result = 0;
      result = mPFunction[fNumber]->Calculate(y) - resultCoefficient;
      return result;
    }
    else
    {
      return mPFunction[fNumber]->Calculate(y);
    }
  }
}


template <class FType>
class TConstrainedProblem : virtual public TConstrainedProblemBase <FType>
{
public:
  TConstrainedProblem()
  {
  };

  TConstrainedProblem(const TConstrainedProblem& problem) : TConstrainedProblemBase<FType>(problem)
  {
  };

  TConstrainedProblem(FType** PFunction, double** q, double* ZoomRatios, double** Shift,
    bool IsImprovementOfTheObjective, double** ImprovementCoefficients, int CountConstrained, int NumberOfCriterions, int Dim) :
  TConstrainedProblemBase<FType>(PFunction, q, ZoomRatios, Shift,
    IsImprovementOfTheObjective, ImprovementCoefficients,
    CountConstrained, NumberOfCriterions, Dim)
  {
  }

  /** ����� ���������� ���������� ����� ����������� �������� ������� �������
  \param[out] y �����, � ������� ����������� ����������� ��������
  \return ��� ������ (#OK ��� #UNDEFINED)
  */
  virtual int GetOptimumPoint(double* point) const
  {
    this->mPFunction[this->mNumberOfConstraints]->GetOptimumPoint(point);

    for (int i = 0; i < this->mDim; i++)
    {
      point[i] = point[i] * this->mZoomRatios[this->mNumberOfConstraints]
        + this->mShift[this->mNumberOfConstraints][i];
    }

    return TProblemPar::ProblemOK;
  }

  /// ���������� ����� ����������� �������� ��� ������� fNumber
  virtual int GetConstraintOptimumPoint(double* point, int fNumber)
  {
    this->mPFunction[fNumber]->GetOptimumPoint(point);

    for (int i = 0; i < this->mDim; i++)
    {
      point[i] = point[i] * this->mZoomRatios[fNumber] + this->mShift[fNumber][i];
    }

    return TProblemPar::ProblemOK;
  }
};


#endif // CONSTRAINED_PROBLEM_H
