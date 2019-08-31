/*
* author: pedro.borges.melo@gmail.com
* date: August/2019
*
* Efficient implementation of the upper smoothing for a two stage 
* stochastic problem with quadratic problems on the second stage.
*/

#ifndef DEFINITION_H_
#define DEFINITION_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#include "IpStdCInterface.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

int solveNonDiffProbUpperSmoothing(
	int maxIterPerSubProblem,
	int totalNewtonRefinementStepsAfterIpopt,
	bool useSecondOrderInfoMasterProblem,
	double epsilon,
	double tikhonovMu,
	int totalDecisionVarsFirstStage,
	int totalScenarios,
	int* totalDecisionVarsPerScenario,
	int totalEqConsFirstStage,
	double* rhsFirstStage,
	int* totalEqConstrPerScenario,
	double* costFirstStage,
	double* probVecScenarios,
	double* costPerScenario,
	double* rhsPerScenario,
	int totalNonZeroFirstStageW,
	int* nzRowsFirstStageW,
	int* nzColsFirstStageW,
	double* nzValsFirstStageW,
	int* totalNonZeroPerScenarioW,
	int* nzRowsPerScenarioW,
	int* nzColsPerScenarioW,
	double* nzValsPerScenarioW,
	int* nzRowsPerScenarioWt,
	int* nzColsPerScenarioWt,
	double* nzValsPerScenarioWt,
	int* totalNonZerPerScenarioT,
	int* nzRowsPerScenarioT,
	int* nzColsPerScenarioT,
	double* nzValsPerScenarioT,
	int* nzRowsPerScenarioTt,
	int* nzColsPerScenarioTt,
	double* nzValsPerScenarioTt,
	int totalNzQuadPerturbation,
	int* nzRowsQuadPerturbation,
	int* nzColsQuadPerturbation,
	double* nzValsQuadPerturbation,	
	int* totalNzQuadPerturbationPerScenario,
	int* nzRowsQuadPerturbationPerScenario,
	int* nzColsQuadPerturbationPerScenario,
	double* nzValsQuadPerturbationPerScenario,
	double* optSolOut,
	double* optValOut,
	double* optMulEqConsOut);

#ifdef __cplusplus
}
#endif

struct QpProbData {
	bool useSecondOrderInfo;
	int maxIterIpopt;
	int totalNewtonRefinementStepsAfterIpopt;
	double muTarget;
	double tikhonovMu;
	int totalDecisionVars;
	int totalEqualityConstraints;
	double* costVec;
	int totalNzQuadCost;
	int* idxRowNzQuadCost;
	int* idxColNzQuadCost;
	double* valNzQuadCost;
	int totalNzSymQuadCost;
	int* idxRowNzSymQuadCost;
	int* idxColNzSymQuadCost;
	double* valNzSymQuadCost;
	int totalNzMatrixA;
	int* idxRowNzMatrixA;
	int* idxColNzMatrixA;
	double* valNzMatrixA;
	int* idxRowNzMatrixAt;
	int* idxColNzMatrixAt;
	double* valNzMatrixAt;
	double* rhsVec;
};

struct pardisoData {
	void *pt[64];
	MKL_INT maxfct;
	MKL_INT mnum;
	MKL_INT mtype;
	MKL_INT phase;
	MKL_INT idum;
	MKL_INT nrhs;
	MKL_INT iparm[64];
	MKL_INT msglvl;
	MKL_INT error;
	double ddum;
};

struct DataSolverUpperSmoothing {
	// control var
	int levelIter;
	int scenarioIter;
	bool useSecondOrderInfoSecondStage;
	// gross data
	double epsilon;
	double tikhonovMu;
	int totalDecisionVarsFirstStage;
	int totalScenarios;
	int* totalDecisionVarsPerScenario;
	int totalEqConsFirstStage;
	double* rhsFirstStage;
	int* totalEqConstrPerScenario;
	double* costFirstStage;
	double* probVecScenarios;
	double* costPerScenario;
	double* rhsPerScenario;
	int totalNonZeroFirstStageW;
	int* nzRowsFirstStageW;
	int* nzColsFirstStageW;
	double* nzValsFirstStageW;
	int* totalNonZeroPerScenarioW;
	int* nzRowsPerScenarioW;
	int* nzColsPerScenarioW;
	double* nzValsPerScenarioW;
	int* nzRowsPerScenarioWt;
	int* nzColsPerScenarioWt;
	double* nzValsPerScenarioWt;
	int* totalNonZerPerScenarioT;
	int* nzRowsPerScenarioT;
	int* nzColsPerScenarioT;
	double* nzValsPerScenarioT;
	int* nzRowsPerScenarioTt;
	int* nzColsPerScenarioTt;
	double* nzValsPerScenarioTt;
	int totalNzQuadPerturbation;
	int* nzRowsQuadPerturbation;
	int* nzColsQuadPerturbation;
	double* nzValsQuadPerturbation;
	int* totalNzQuadPerturbationPerScenario;
	int* nzRowsQuadPerturbationPerScenario;
	int* nzColsQuadPerturbationPerScenario;
	double* nzValsQuadPerturbationPerScenario;
	// organized data
	struct QpProbData firstStageData;
	std::vector<struct QpProbData> secondStageData;
	// pre-computed newton matrices per scenario
	// NOTE: stored in symmetric form
	// NOTE: has to take care with the sing of the multiplier
	std::vector<int> sizesNewtonSystemsPerScenario;
	std::vector<int*> rowIncrementsNewtonSystemsPerScenario;
	std::vector<int*> colsNewtonSystemsPerScenario;
	std::vector<double*> valsNewtonSystemsPerScenario;
	std::vector<int*> positionsOfDiagonalElementsToChangePerScenario;
	double* rhsForNewtonPerScenario;
	double* solForNewtonPerScenario;
	struct pardisoData pardisoData_;
	// buffer for solutions of scenario problems	
	double* optPrimalSolScen;
	double* optMultEqConstraintsScen;
	// buffer derivatives of solution mappings of scenario problems
	// NOTE: the size of this vectors total decision vars of the first stage
	// NOTE: for our purposes there is no need to store the dual solutions
	std::vector<double*> derivativesPrimalPerScenarioWrtPcoord;
	// buffer for agregated values of second stage problem
	double* lastFirstStageTrial;
	double lastValuesPerScenario;
	double* lastGradSecondStage;
	int totalPositionsHessSecondStage;
	int* idxRowHessSecondStage;
	int* idxColHessSecondStage;
	double* valsHessSecondStage;
	// buffer for refinement of solution of scenario problem
	double* gradientOfLagrangianIter;
	double* gradLagrTrialPoint;
	double* newtonDirectionIter;
	double* primalSolIter;
	double* primalSolIterDamp;
	double* multiplierIter; 
	double* multIterDamp;
};

double normVec(double* vec, int size);

int initPardisoInternalDataForSymmetricLinearSystem(
	struct pardisoData* pardisoData_,
	int n, int* ia, int* ja, double* a);

int solveSymmetricSystemWithPardiso(
	int n, int* ia, int* ja, double* a,
	double* b, double* x);

int destroyPardisoForSymmetricSystem(
	struct pardisoData* pardisoData_,
	int n, int* ia, int* ja);

int allSolveSymmetricSystemWithPardiso(
	int n, int* ia, int* ja, double* a,
	double* b, double* x);

int getStartingPointForIpoptWithCplex(
	struct QpProbData* data,
	double* initialPrimalIterate,
	double* initialDualIterate);

void computeGradientOfLagrangian(
	double* gradientOfLagrangian,
	struct QpProbData* data,
	double* primalSolIter,
	double* multiplierIter);

int refineMuTargetSolutionWithDampedNewtonMethod(
	int scenario,
	struct DataSolverUpperSmoothing* user_data,
	double* initialPrimalIterate,
	double* initialDualIterate);

int solveMuTargetSubproblemWithIpopt(
	int scenario,
	struct DataSolverUpperSmoothing* user_data,
	double* initialPrimalIterate,
	double* optVal,
	double* initialDualIterate);

void updateSecondStageInformationIfNeeded(
	struct DataSolverUpperSmoothing* user_data,
	Number* x_trial);

void shiftVectorToPositiveOrthantIfNeeded(
	int size, double* vec, double MIN_VALUE);

Bool eval_f_quad_prob(Index n, Number* x, Bool new_x, Number* obj_value, UserDataPtr user_data_all);

Bool eval_grad_f_quad_prob(Index n, Number* x, Bool new_x, Number* grad_f, UserDataPtr user_data_all);

Bool eval_g_quad_prob(Index n, Number* x, Bool new_x, Index m, Number* g, UserDataPtr user_data_all);

Bool eval_jac_g_quad_prob(Index n, Number *x, Bool new_x, Index m, Index nele_jac, 
	Index *iRow, Index *jCol, Number *values, UserDataPtr user_data_all);

Bool eval_h_quad_prob(Index n, Number *x, Bool new_x, Number obj_factor, Index m, 
	Number *lambda, Bool new_lambda, Index nele_hess, 
	Index *iRow, Index *jCol, Number *values, UserDataPtr user_data_all);

void initDataSolverUpperSmoothing(
	DataSolverUpperSmoothing* data_, 
	bool useSecondOrderInfoMasterProblem,
	int maxIterPerSubProblem,
	int totalNewtonRefinementStepsAfterIpopt_,
	double epsilon_,
	double tikhonovMu_,
	int totalDecisionVarsFirstStage_,
	int totalScenarios_,
	int* totalDecisionVarsPerScenario_,
	int totalEqConsFirstStage_,
	double* rhsFirstStage_,
	int* totalEqConstrPerScenario_,
	double* costFirstStage_,
	double* probVecScenarios_,
	double* costPerScenario_,
	double* rhsPerScenario_,
	int totalNonZeroFirstStageW_,
	int* nzRowsFirstStageW_,
	int* nzColsFirstStageW_,
	double* nzValsFirstStageW_,
	int* totalNonZeroPerScenarioW_,
	int* nzRowsPerScenarioW_,
	int* nzColsPerScenarioW_,
	double* nzValsPerScenarioW_,
	int* nzRowsPerScenarioWt_,
	int* nzColsPerScenarioWt_,
	double* nzValsPerScenarioWt_,
	int* totalNonZerPerScenarioT_,
	int* nzRowsPerScenarioT_,
	int* nzColsPerScenarioT_,
	double* nzValsPerScenarioT_,
	int* nzRowsPerScenarioTt_,
	int* nzColsPerScenarioTt_,
	double* nzValsPerScenarioTt_,
	int totalNzQuadPerturbation_,
	int* nzRowsQuadPerturbation_,
	int* nzColsQuadPerturbation_,
	double* nzValsQuadPerturbation_,
	int* totalNzQuadPerturbationPerScenario_,
	int* nzRowsQuadPerturbationPerScenario_,
	int* nzColsQuadPerturbationPerScenario_,
	double* nzValsQuadPerturbationPerScenario_);

void freeDataSolverUpperSmoothing(
	DataSolverUpperSmoothing* data_);

#endif
