/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2016 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      problemPar.h                                         //
//                                                                         //
//  Purpose:   Header file for ExaMin problem interface                    //
//                                                                         //
//                                                                         //
//  Author(s): Lebedev I.                                                 //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/**
\file problemPar.h

\authors Лебедев И.
\date 2016
\copyright ННГУ им. Н.И. Лобачевского

\brief Объявление абстрактного класса #TProblemPar

\details Объявление абстрактного класса #TProblemPar, являющегося одним из базовых для параметризованных задач
*/

#ifndef __PROBLEM_PAR_H__
#define __PROBLEM_PAR_H__

class TProblemPar
{
protected:
  /// Число критериев
  int mNumberOfCriterions;
  /// Число ограничений
  int mNumberOfConstraints;
  /// Разерность задачи
  int mDim;


  /** Метод возвращает число общее функций в задаче (оно равно число ограничений + число критериев)
  \return Число функций
  */
  virtual int GetRealNumberOfFunctions() const
  {
    return mNumberOfCriterions + mNumberOfConstraints;
  }

  /** Метод возвращает число ограничений в задаче
  \return Число ограничений
  */
  virtual int GetRealNumberOfConstraints() const
  {
    return mNumberOfConstraints;
  }
public:


#ifndef AccuracyDouble
#define AccuracyDouble 0.00000001
#endif
#ifndef MAX_DIM
#define MAX_DIM 50
#endif

  /// Код ошибки, возвращаемый, если операция завершена успешно
  static const int ProblemOK = 0;
  /** Код ошибки, возвращаемый методами #GetOptimumValue и #GetOptimumPoint,
  если соответствующие параметры задачи не определены,
  */
  static const int ProblemUNDEFINED = -1;
  /// Код ошибки, возвращаемый, если операция не выполнена
  static const int ProblemERROR = -2;

  /// Возвращает точку глобального оптимума для функции fNumber
  virtual int GetConstraintOptimumPoint(double* point, int fNumber) = 0;

  /** Метод возвращает координаты точки глобального минимума целевой функции
  \param[out] y точка, в которой достигается оптимальное значение
  \return Код ошибки (#OK или #UNDEFINED)
  */
  virtual int GetOptimumPoint(double* y) const  = 0;

  TProblemPar()
  {
    mNumberOfCriterions = 0;
    mNumberOfConstraints = 0;
    mDim = 0;
  }
};


#endif
// - end of file ----------------------------------------------------------------------------------
