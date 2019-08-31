#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "2stpUpperSmoothing.h"

int test01() {
	printf(" \n\n >>>>>>>>>> test 01 >>>>>>>>>> \n\n ");
	// data
	bool useSecondOrderInfoMasterProblem = true;
	double epsilon = 0.5;
	double tikhonovMu = 0.0;
	int totalDecisionVarsFirstStage = 3;
	int totalScenarios = 2;
	int totalDecisionVarsPerScenario[] = {3, 3};
	int totalEqConsFirstStage = 2;
	double rhsFirstStage[] = {1.0, 2.0};
	int totalEqConstrPerScenario[] = {2, 2};
	double costFirstStage[] = {1.0, 2.0, 1.0};
	double probVecScenarios[] = {0.5, 0.5};
	double costPerScenario[] = {1.0, 2.0, 1.0, 1.0, 2.0, 1.0};
	double rhsPerScenario[] = {1, 2, 1, 2};
	int totalNonZeroFirstStageW = 2;
	int nzRowsFirstStageW[] = {1, 2};
	int nzColsFirstStageW[] = {1, 2};
	double nzValsFirstStageW[] = {1, 1};
	int totalNonZeroPerScenarioW[] = {2, 2};
	int nzRowsPerScenarioW[] = {1, 2, 1, 2};
	int nzColsPerScenarioW[] = {1, 2, 1, 2};
	double nzValsPerScenarioW[] = {1, 1, 1, 1};
	int nzRowsPerScenarioWt[] = {1, 2, 1, 2};
	int nzColsPerScenarioWt[] = {1, 2, 1, 2};
	double nzValsPerScenarioWt[] = {1, 1, 1, 1};
	int totalNonZerPerScenarioT[] = {0, 0};
	int* nzRowsPerScenarioT = NULL;
	int* nzColsPerScenarioT = NULL;
	double* nzValsPerScenarioT = NULL;
	int* nzRowsPerScenarioTt = NULL;
	int* nzColsPerScenarioTt = NULL;
	double* nzValsPerScenarioTt = NULL;
	int totalNzQuadPerturbation = 0;
	int* nzRowsQuadPerturbation = NULL;
	int* nzColsQuadPerturbation = NULL;
	double* nzValsQuadPerturbation = NULL;
	int totalNzQuadPerturbationPerScenario[] = {3, 3};
	int nzRowsQuadPerturbationPerScenario[] = {1, 2, 3, 1, 2, 3};
	int nzColsQuadPerturbationPerScenario[] = {1, 2, 3, 1, 2, 3};
	double nzValsQuadPerturbationPerScenario[] = {0, 0, 0, 0, 0, 0};
	double optValOut = 5.0;
	double optSolOut[] = {1.0, 2.0, 0.0};
	double optMulEqConsOut[] = {0.0, 0.0};
	int maxIterPerSubProblem = 100;
	int totalNewtonRefinementStepsAfterIpopt = 4;
	// solve
	solveNonDiffProbUpperSmoothing(
		maxIterPerSubProblem,
		totalNewtonRefinementStepsAfterIpopt,
		useSecondOrderInfoMasterProblem,
		epsilon,
		tikhonovMu,
		totalDecisionVarsFirstStage,
		totalScenarios,
		totalDecisionVarsPerScenario,
		totalEqConsFirstStage,
		rhsFirstStage,
		totalEqConstrPerScenario,
		costFirstStage,
		probVecScenarios,
		costPerScenario,
		rhsPerScenario,
		totalNonZeroFirstStageW,
		nzRowsFirstStageW,
		nzColsFirstStageW,
		nzValsFirstStageW,
		totalNonZeroPerScenarioW,
		nzRowsPerScenarioW,
		nzColsPerScenarioW,
		nzValsPerScenarioW,
		nzRowsPerScenarioWt,
		nzColsPerScenarioWt,
		nzValsPerScenarioWt,
		totalNonZerPerScenarioT,
		nzRowsPerScenarioT,
		nzColsPerScenarioT,
		nzValsPerScenarioT,
		nzRowsPerScenarioTt,
		nzColsPerScenarioTt,
		nzValsPerScenarioTt,
		totalNzQuadPerturbation,
		nzRowsQuadPerturbation,
		nzColsQuadPerturbation,
		nzValsQuadPerturbation,
		totalNzQuadPerturbationPerScenario,
		nzRowsQuadPerturbationPerScenario,
		nzColsQuadPerturbationPerScenario,
		nzValsQuadPerturbationPerScenario,
		optSolOut,
		&optValOut,
		optMulEqConsOut);
	// print
	printf("\n");
	printf(" optValOut: %f \n", optValOut);
	printf(" correct optValOut = 10 \n");
	printf("\n");
	for (int i=0; i<totalDecisionVarsFirstStage; ++i) {
		printf(" optSolOut[%d] = %f \n", i, optSolOut[i]);
	}
	printf("\n");
	for (int i=0; i<totalEqConsFirstStage; ++i) {
		printf(" optMulEqConsOut[%d] = %f \n", i, optMulEqConsOut[i]);
	}
	// return
	return 0;
}

/*
* min x + Q(x) s.t. x = 0.555, x >=0 AND Q(x) = min { y : y>=0, y = 1 - x}
*/

int test02() {
	printf(" \n\n >>>>>>>>>> test 02 >>>>>>>>>> \n\n ");
	// data
	bool useSecondOrderInfoMasterProblem = true;
	double epsilon = 0.001;
	double tikhonovMu = 0.0;
	int totalNzQuadPerturbation = 0;
	int* nzRowsQuadPerturbation = NULL;
	int* nzColsQuadPerturbation = NULL;
	double* nzValsQuadPerturbation = NULL;
	int totalDecisionVarsFirstStage = 1;
	int totalScenarios = 1;
	int totalDecisionVarsPerScenario[] = {1};
	int totalEqConsFirstStage = 1;
	double rhsFirstStage[] = {0.555};
	int totalEqConstrPerScenario[] = {1};
	double costFirstStage[] = {1.0};
	double probVecScenarios[] = {1.0};
	double costPerScenario[] =  {1.0};
	double rhsPerScenario[] = {1};
	int totalNonZeroFirstStageW = 1;
	int nzRowsFirstStageW[] = {1};
	int nzColsFirstStageW[] = {1};
	double nzValsFirstStageW[] = {1};
	int totalNonZeroPerScenarioW[] = {1};
	int nzRowsPerScenarioW[] = {1};
	int nzColsPerScenarioW[] = {1};
	double nzValsPerScenarioW[] = {1};
	int nzRowsPerScenarioWt[] = {1};
	int nzColsPerScenarioWt[] = {1};
	double nzValsPerScenarioWt[] = {1};
	int totalNonZerPerScenarioT[] = {1};
	int nzRowsPerScenarioT[] = {1};
	int nzColsPerScenarioT[] = {1};
	double nzValsPerScenarioT[] = {1};
	int nzRowsPerScenarioTt[] = {1};
	int nzColsPerScenarioTt[] = {1};
	double nzValsPerScenarioTt[] = {1};
	int totalNzQuadPerturbationPerScenario[] = {1};
	int nzRowsQuadPerturbationPerScenario[] = {1};
	int nzColsQuadPerturbationPerScenario[] = {1};
	double nzValsQuadPerturbationPerScenario[] = {0.0};
	double optSolOut[] = {0.0};
	double optValOut = 0.0;
	double optMulEqConsOut[] = {0.0};
	int maxIterPerSubProblem = 100;
	int totalNewtonRefinementStepsAfterIpopt = 4;
	// solve
	solveNonDiffProbUpperSmoothing(
		maxIterPerSubProblem,
		totalNewtonRefinementStepsAfterIpopt,
		useSecondOrderInfoMasterProblem,
		epsilon,
		tikhonovMu,
		totalDecisionVarsFirstStage,
		totalScenarios,
		totalDecisionVarsPerScenario,
		totalEqConsFirstStage,
		rhsFirstStage,
		totalEqConstrPerScenario,
		costFirstStage,
		probVecScenarios,
		costPerScenario,
		rhsPerScenario,
		totalNonZeroFirstStageW,
		nzRowsFirstStageW,
		nzColsFirstStageW,
		nzValsFirstStageW,
		totalNonZeroPerScenarioW,
		nzRowsPerScenarioW,
		nzColsPerScenarioW,
		nzValsPerScenarioW,
		nzRowsPerScenarioWt,
		nzColsPerScenarioWt,
		nzValsPerScenarioWt,
		totalNonZerPerScenarioT,
		nzRowsPerScenarioT,
		nzColsPerScenarioT,
		nzValsPerScenarioT,
		nzRowsPerScenarioTt,
		nzColsPerScenarioTt,
		nzValsPerScenarioTt,
		totalNzQuadPerturbation,
		nzRowsQuadPerturbation,
		nzColsQuadPerturbation,
		nzValsQuadPerturbation,
		totalNzQuadPerturbationPerScenario,
		nzRowsQuadPerturbationPerScenario,
		nzColsQuadPerturbationPerScenario,
		nzValsQuadPerturbationPerScenario,
		optSolOut,
		&optValOut,
		optMulEqConsOut);
	// print
	printf("\n");
	printf(" optValOut: %f \n", optValOut);
	printf(" correct optValOut = 1\n", optValOut);
	printf("\n");
	for (int i=0; i<totalDecisionVarsFirstStage; ++i) {
		printf(" optSolOut[%d] = %f \n", i, optSolOut[i]);
	}
	printf("\n");
	for (int i=0; i<totalEqConsFirstStage; ++i) {
		printf(" optMulEqConsOut[%d] = %f \n", i, optMulEqConsOut[i]);
	}
	// return
	return 0;
}

/*
* min {Q(x) : x >= 0} AND Q(x) = min {y : y >= 0, y >= |5-x|}
*/

int test05() {
	printf(" \n\n >>>>>>>>>> test 05 >>>>>>>>>> \n\n ");
	// first stage data
	bool useSecondOrderInfoMasterProblem = true;
	double epsilon = 0.000001;
	double tikhonovMu = 0.0;
	int totalNzQuadPerturbation = 0;
	int* nzRowsQuadPerturbation = NULL;
	int* nzColsQuadPerturbation = NULL;
	double* nzValsQuadPerturbation = NULL;
	int totalDecisionVarsFirstStage = 1;
	int totalEqConsFirstStage = 0;
	double* rhsFirstStage = NULL;
	double costFirstStage[] = {0.0};
	int totalNonZeroFirstStageW = 0;
	int* nzRowsFirstStageW = NULL;
	int* nzColsFirstStageW = NULL;
	double* nzValsFirstStageW = NULL;
	// second stage data
	int totalScenarios = 1;
	int totalDecisionVarsPerScenario[] = {3};
	int totalEqConstrPerScenario[] = {2};
	double probVecScenarios[] = {1.0};
	double costPerScenario[] =  {1.0, 0.0, 0.0};
	double rhsPerScenario[] = {5.0, -5.0};
	int totalNonZeroPerScenarioW[] = {4};
	int nzRowsPerScenarioW[] = {1, 1, 2, 2};
	int nzColsPerScenarioW[] = {1, 2, 1, 3};
	double nzValsPerScenarioW[] = {1.0, -1.0, 1.0, -1.0};
	int nzRowsPerScenarioWt[] = {1, 1, 2, 3};
	int nzColsPerScenarioWt[] = {1, 2, 1, 2};
	double nzValsPerScenarioWt[] = {1.0, 1.0, -1.0, -1.0};
	int totalNonZerPerScenarioT[] = {2};
	int nzRowsPerScenarioT[] = {1, 2};
	int nzColsPerScenarioT[] = {1, 1};
	double nzValsPerScenarioT[] = {1.0, -1.0};
	int nzRowsPerScenarioTt[] = {1, 1};
	int nzColsPerScenarioTt[] = {1, 2};
	double nzValsPerScenarioTt[] = {1.0, -1.0};
	int totalNzQuadPerturbationPerScenario[] = {3};
	int nzRowsQuadPerturbationPerScenario[] = {1, 2, 3};
	int nzColsQuadPerturbationPerScenario[] = {1, 2, 3};
	double nzValsQuadPerturbationPerScenario[] = {0.0, 0.0, 0.0};
	double optValOut = 0.0;
	double optSolOut[] = {0.0};
	double* optMulEqConsOut = NULL;
	int maxIterPerSubProblem = 100;
	int totalNewtonRefinementStepsAfterIpopt = 4;
	// solve
	solveNonDiffProbUpperSmoothing(
		maxIterPerSubProblem,
		totalNewtonRefinementStepsAfterIpopt,
		useSecondOrderInfoMasterProblem,
		epsilon,
		tikhonovMu,
		totalDecisionVarsFirstStage,
		totalScenarios,
		totalDecisionVarsPerScenario,
		totalEqConsFirstStage,
		rhsFirstStage,
		totalEqConstrPerScenario,
		costFirstStage,
		probVecScenarios,
		costPerScenario,
		rhsPerScenario,
		totalNonZeroFirstStageW,
		nzRowsFirstStageW,
		nzColsFirstStageW,
		nzValsFirstStageW,
		totalNonZeroPerScenarioW,
		nzRowsPerScenarioW,
		nzColsPerScenarioW,
		nzValsPerScenarioW,
		nzRowsPerScenarioWt,
		nzColsPerScenarioWt,
		nzValsPerScenarioWt,
		totalNonZerPerScenarioT,
		nzRowsPerScenarioT,
		nzColsPerScenarioT,
		nzValsPerScenarioT,
		nzRowsPerScenarioTt,
		nzColsPerScenarioTt,
		nzValsPerScenarioTt,
		totalNzQuadPerturbation,
		nzRowsQuadPerturbation,
		nzColsQuadPerturbation,
		nzValsQuadPerturbation,
		totalNzQuadPerturbationPerScenario,
		nzRowsQuadPerturbationPerScenario,
		nzColsQuadPerturbationPerScenario,
		nzValsQuadPerturbationPerScenario,
		optSolOut,
		&optValOut,
		optMulEqConsOut);
	// print
	printf("\n");
	printf(" optValOut: %f \n", optValOut);
	printf(" correct optValOut = 0.0 \n");
	printf("\n");
	for (int i=0; i<totalDecisionVarsFirstStage; ++i) {
		printf(" optSolOut[%d] = %f \n", i, optSolOut[i]);
	}
	printf("\n");
	for (int i=0; i<totalEqConsFirstStage; ++i) {
		printf(" optMulEqConsOut[%d] = %f \n", i, optMulEqConsOut[i]);
	}
	// return
	return 0;
}

/*
* min {Q(x) : x >= 0} AND Q(x) = min {y + y^2 : y >= 0, y >= |5-x|}
*/

int test06() {
	printf(" \n\n >>>>>>>>>> test 06 >>>>>>>>>> \n\n ");
	// first stage data
	bool useSecondOrderInfoMasterProblem = true;
	double epsilon = 0.001;
	double tikhonovMu = 0.0;
	int totalNzQuadPerturbation = 0;
	int* nzRowsQuadPerturbation = NULL;
	int* nzColsQuadPerturbation = NULL;
	double* nzValsQuadPerturbation = NULL;
	int totalDecisionVarsFirstStage = 1;
	int totalEqConsFirstStage = 0;
	double* rhsFirstStage = NULL;
	double costFirstStage[] = {0.0};
	int totalNonZeroFirstStageW = 0;
	int* nzRowsFirstStageW = NULL;
	int* nzColsFirstStageW = NULL;
	double* nzValsFirstStageW = NULL;
	// second stage data
	int totalScenarios = 1;
	int totalDecisionVarsPerScenario[] = {3};
	int totalEqConstrPerScenario[] = {2};
	double probVecScenarios[] = {1.0};
	double costPerScenario[] =  {1.0, 0.0, 0.0};
	double rhsPerScenario[] = {5.0, -5.0};
	int totalNonZeroPerScenarioW[] = {4};
	int nzRowsPerScenarioW[] = {1, 1, 2, 2};
	int nzColsPerScenarioW[] = {1, 2, 1, 3};
	double nzValsPerScenarioW[] = {1.0, -1.0, 1.0, -1.0};
	int nzRowsPerScenarioWt[] = {1, 1, 2, 3};
	int nzColsPerScenarioWt[] = {1, 2, 1, 2};
	double nzValsPerScenarioWt[] = {1.0, 1.0, -1.0, -1.0};
	int totalNonZerPerScenarioT[] = {2};
	int nzRowsPerScenarioT[] = {1, 2};
	int nzColsPerScenarioT[] = {1, 1};
	double nzValsPerScenarioT[] = {1.0, -1.0};
	int nzRowsPerScenarioTt[] = {1, 1};
	int nzColsPerScenarioTt[] = {1, 2};
	double nzValsPerScenarioTt[] = {1.0, -1.0};
	int totalNzQuadPerturbationPerScenario[] = {3};
	int nzRowsQuadPerturbationPerScenario[] = {1, 2, 3};
	int nzColsQuadPerturbationPerScenario[] = {1, 2, 3};
	double nzValsQuadPerturbationPerScenario[] = {1.0, 0.0, 0.0};
	double optValOut = 0.0;
	double optSolOut[] = {0.0};
	double* optMulEqConsOut = NULL;
	int maxIterPerSubProblem = 100;
	int totalNewtonRefinementStepsAfterIpopt = 4;
	// solve
	solveNonDiffProbUpperSmoothing(
		maxIterPerSubProblem,
		totalNewtonRefinementStepsAfterIpopt,
		useSecondOrderInfoMasterProblem,
		epsilon,
		tikhonovMu,
		totalDecisionVarsFirstStage,
		totalScenarios,
		totalDecisionVarsPerScenario,
		totalEqConsFirstStage,
		rhsFirstStage,
		totalEqConstrPerScenario,
		costFirstStage,
		probVecScenarios,
		costPerScenario,
		rhsPerScenario,
		totalNonZeroFirstStageW,
		nzRowsFirstStageW,
		nzColsFirstStageW,
		nzValsFirstStageW,
		totalNonZeroPerScenarioW,
		nzRowsPerScenarioW,
		nzColsPerScenarioW,
		nzValsPerScenarioW,
		nzRowsPerScenarioWt,
		nzColsPerScenarioWt,
		nzValsPerScenarioWt,
		totalNonZerPerScenarioT,
		nzRowsPerScenarioT,
		nzColsPerScenarioT,
		nzValsPerScenarioT,
		nzRowsPerScenarioTt,
		nzColsPerScenarioTt,
		nzValsPerScenarioTt,
		totalNzQuadPerturbation,
		nzRowsQuadPerturbation,
		nzColsQuadPerturbation,
		nzValsQuadPerturbation,
		totalNzQuadPerturbationPerScenario,
		nzRowsQuadPerturbationPerScenario,
		nzColsQuadPerturbationPerScenario,
		nzValsQuadPerturbationPerScenario,
		optSolOut,
		&optValOut,
		optMulEqConsOut);
	// print
	printf("\n");
	printf(" optValOut: %f \n", optValOut);
	printf(" correct optValOut = 0.0 \n");
	printf("\n");
	for (int i=0; i<totalDecisionVarsFirstStage; ++i) {
		printf(" optSolOut[%d] = %f \n", i, optSolOut[i]);
	}
	printf(" correct optSolOut[0] = 5.0 \n");
	printf("\n");
	for (int i=0; i<totalEqConsFirstStage; ++i) {
		printf(" optMulEqConsOut[%d] = %f \n", i, optMulEqConsOut[i]);
	}
	// return
	return 0;
}

/*
* min {Q(x_1, x_2) : x_1, x_2 >= 0} AND Q(x_1, x_2) = min {y_1 + y_2: y_1, y_2 >= 0, y_1 >= |5-x_1|, y_2 >= |5-x_2|}
*/
int test07() {
	printf(" \n\n >>>>>>>>>> test 07 >>>>>>>>>> \n\n ");
	// first stage data
	bool useSecondOrderInfoMasterProblem = true;
	double epsilon = 0.0001;
	double tikhonovMu = 1.0;
	int totalNzQuadPerturbation = 0;
	int* nzRowsQuadPerturbation = NULL;
	int* nzColsQuadPerturbation = NULL;
	double* nzValsQuadPerturbation = NULL;
	int totalDecisionVarsFirstStage = 2;
	int totalEqConsFirstStage = 0;
	double* rhsFirstStage = NULL;
	double costFirstStage[] = {0.0, 0.0};
	int totalNonZeroFirstStageW = 0;
	int* nzRowsFirstStageW = NULL;
	int* nzColsFirstStageW = NULL;
	double* nzValsFirstStageW = NULL;
	// second stage data
	int totalScenarios = 1;
	int totalDecisionVarsPerScenario[] = {6};
	int totalEqConstrPerScenario[] = {4};
	double probVecScenarios[] = {1.0};
	double costPerScenario[] =  {1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
	double rhsPerScenario[] = {5.0, -5.0, 5.0, -5.0};
	int totalNonZeroPerScenarioW[] = {8};
	int nzRowsPerScenarioW[] = {1, 1, 2, 2, 3, 3, 4, 4};
	int nzColsPerScenarioW[] = {1, 2, 1, 3, 4, 5, 4, 6};
	double nzValsPerScenarioW[] = {1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
	int nzRowsPerScenarioWt[] = {1, 1, 2, 3, 4, 4, 5, 6};
	int nzColsPerScenarioWt[] = {1, 2, 1, 2, 3, 4, 3, 4};
	double nzValsPerScenarioWt[] = {1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0};
	int totalNonZerPerScenarioT[] = {4};
	int nzRowsPerScenarioT[] = {1, 2, 3, 4};
	int nzColsPerScenarioT[] = {1, 1, 2, 2};
	double nzValsPerScenarioT[] = {1.0, -1.0, 1.0, -1.0};
	int nzRowsPerScenarioTt[] = {1, 1, 2, 2};
	int nzColsPerScenarioTt[] = {1, 2, 3, 4};
	double nzValsPerScenarioTt[] = {1.0, -1.0, 1.0, -1.0};
	int totalNzQuadPerturbationPerScenario[] = {0};
	int* nzRowsQuadPerturbationPerScenario = NULL;
	int* nzColsQuadPerturbationPerScenario = NULL;
	double* nzValsQuadPerturbationPerScenario = NULL;
	double optValOut = 0.0;
	double optSolOut[] = {0.0, 0.0};
	double* optMulEqConsOut = NULL;
	int maxIterPerSubProblem = 100;
	int totalNewtonRefinementStepsAfterIpopt = 0;
	// solve
	solveNonDiffProbUpperSmoothing(
		maxIterPerSubProblem,
		totalNewtonRefinementStepsAfterIpopt,
		useSecondOrderInfoMasterProblem,
		epsilon,
		tikhonovMu,
		totalDecisionVarsFirstStage,
		totalScenarios,
		totalDecisionVarsPerScenario,
		totalEqConsFirstStage,
		rhsFirstStage,
		totalEqConstrPerScenario,
		costFirstStage,
		probVecScenarios,
		costPerScenario,
		rhsPerScenario,
		totalNonZeroFirstStageW,
		nzRowsFirstStageW,
		nzColsFirstStageW,
		nzValsFirstStageW,
		totalNonZeroPerScenarioW,
		nzRowsPerScenarioW,
		nzColsPerScenarioW,
		nzValsPerScenarioW,
		nzRowsPerScenarioWt,
		nzColsPerScenarioWt,
		nzValsPerScenarioWt,
		totalNonZerPerScenarioT,
		nzRowsPerScenarioT,
		nzColsPerScenarioT,
		nzValsPerScenarioT,
		nzRowsPerScenarioTt,
		nzColsPerScenarioTt,
		nzValsPerScenarioTt,
		totalNzQuadPerturbation,
		nzRowsQuadPerturbation,
		nzColsQuadPerturbation,
		nzValsQuadPerturbation,
		totalNzQuadPerturbationPerScenario,
		nzRowsQuadPerturbationPerScenario,
		nzColsQuadPerturbationPerScenario,
		nzValsQuadPerturbationPerScenario,
		optSolOut,
		&optValOut,
		optMulEqConsOut);
	// print
	printf("\n");
	printf(" optValOut: %f \n", optValOut);
	printf(" correct optValOut = 0.0 \n");
	printf("\n");
	for (int i=0; i<totalDecisionVarsFirstStage; ++i) {
		printf(" optSolOut[%d] = %f \n", i, optSolOut[i]);
	}
	printf(" correct optSolOut[0] = 5.0 \n");
	printf(" correct optSolOut[1] = 5.0 \n");
	printf("\n");
	for (int i=0; i<totalEqConsFirstStage; ++i) {
		printf(" optMulEqConsOut[%d] = %f \n", i, optMulEqConsOut[i]);
	}
	// return
	return 0;
}

int main() {
	// run test
	test01();
	test02();
	test05();
	test06();
	test07();
	// return
	return 0;
}
