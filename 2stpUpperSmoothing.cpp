/*
* author: pedro.borges.melo@gmail.com
* date: August/2019
*
* Efficient implementation of the upper smoothing for a two stage 
* stochastic problem with quadratic problems on the second stage.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#include "IpStdCInterface.h"
#include "2stpUpperSmoothing.h"
#include <ilcp/cp.h>
#include <ilcplex/ilocplexi.h>

using namespace std;

double normVec(double* vec, int size) {
	// init
	double sqrtNorm = 0.0;
	// iteratesolveSmoothedMasterProblem
	for (int i=0; i<size; ++i) {
		sqrtNorm += (vec[i])*(vec[i]);
	}
	// sqrt
	sqrtNorm = sqrt(sqrtNorm);
	// return
	return sqrtNorm;
}

void shiftVectorToPositiveOrthantIfNeeded(int size, double* vec, double MIN_VALUE) {
	for (int i=0; i<size; ++i) {
		if (vec[i] <= MIN_VALUE) {
			vec[i] = MIN_VALUE;
		}
	}
}

int initPardisoInternalDataForSymmetricLinearSystem(
		struct pardisoData* pardisoData_,
		int n, int* ia, int* ja, double* a) {
	// symmetric indefinite matrix
	pardisoData_->mtype = -2;
	// Number of right hand sides
	pardisoData_->nrhs = 1;
	// Setup Pardiso control parameters
	MKL_INT i;
	for ( i = 0; i < 64; i++ ) {
		pardisoData_->iparm[i] = 0;
	}
	// No solver default 
	pardisoData_->iparm[0] = 1;
	// Fill-in reordering from METIS 
	pardisoData_->iparm[1] = 2;
	// No iterative-direct algorithm 
	pardisoData_->iparm[3] = 0;
	// No user fill-in reducing permutation 
	pardisoData_->iparm[4] = 0;
	// Write solution into x 
	pardisoData_->iparm[5] = 0;
	// Not in use 
	pardisoData_->iparm[6] = 0;
	// Max numbers of iterative refinement steps 
	pardisoData_->iparm[7] = 2;
	// Not in use
	pardisoData_->iparm[8] = 0;
	// Perturb the pivot elements with 1E-13 
	pardisoData_->iparm[9] = 13;
	// Use nonsymmetric permutation and scaling MPS 
	pardisoData_->iparm[10] = 1;
	// Not in use 
	pardisoData_->iparm[11] = 0;
	// Maximum weighted matching algorithm is switched-off (default for symmetric)
	// Try iparm[12] = 1 in case of inappropriate accuracy
	pardisoData_->iparm[12] = 0;
	// Output: Number of perturbed pivots 
	pardisoData_->iparm[13] = 0;
	// Not in use 
	pardisoData_->iparm[14] = 0;
	// Not in use 
	pardisoData_->iparm[15] = 0;
	// Not in use 
	pardisoData_->iparm[16] = 0;
	// Output: Number of nonzeros in the factor LU 
	pardisoData_->iparm[17] = -1;
	// Output: Mflops for LU factorization 
	pardisoData_->iparm[18] = -1;
	// Output: Numbers of CG Iterations 
	pardisoData_->iparm[19] = 0;
	// Maximum number of numerical factorizations. 
	pardisoData_->maxfct = 1;
	// Which factorization to use. 
	pardisoData_->mnum = 1;
	// Print statistical information in file 
	pardisoData_->msglvl = 0;
	// Initialize error flag 
	pardisoData_->error = 0;
	//
	// Initialize the internal solver memory pointer. This is only
	// necessary for the FIRST call of the PARDISO solver.
	//
	for ( i = 0; i < 64; i++ ) {
		pardisoData_->pt[i] = 0;
	}
	//
	// Reordering and Symbolic Factorization. This step also allocates
	// all memory that is necessary for the factorization.
	//
	pardisoData_->phase = 11;
	PARDISO (pardisoData_->pt, &(pardisoData_->maxfct), &(pardisoData_->mnum), 
		&(pardisoData_->mtype), &(pardisoData_->phase),
		&n, a, ia, ja, &(pardisoData_->idum), 
		&(pardisoData_->nrhs), pardisoData_->iparm, 
		&(pardisoData_->msglvl), &(pardisoData_->ddum), 
		&(pardisoData_->ddum), &(pardisoData_->error));
	if ( pardisoData_->error != 0 ) {
		printf ("\nERROR during symbolic factorization: %d", pardisoData_->error);
	}
	// printf ("\nReordering completed ... ");
	// printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
	// printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
	//
	// Numerical factorization
	// 
	pardisoData_->phase = 22;
	PARDISO (pardisoData_->pt, &(pardisoData_->maxfct), &(pardisoData_->mnum), 
		&(pardisoData_->mtype), &(pardisoData_->phase),
		&n, a, ia, ja, &(pardisoData_->idum), 
		&(pardisoData_->nrhs), (pardisoData_->iparm), &(pardisoData_->msglvl),
		&(pardisoData_->ddum), &(pardisoData_->ddum), &(pardisoData_->error));
	// Check for error
	if ( pardisoData_->error != 0 ) {
		printf ("\nERROR during numerical factorization: %d", pardisoData_->error);
	}
	// return
	return 0;
}

int solveSymmetricSystemWithPardiso(
		struct pardisoData* pardisoData_,
		int n, int* ia, int* ja, double* a,
		double* b, double* x) {
	// printf ("\nFactorization completed ... ");
	// 
	// Back substitution and iterative refinement
	// 
	pardisoData_->phase = 33;
	// Max numbers of iterative refinement steps
	pardisoData_->iparm[7] = 2;
	PARDISO (pardisoData_->pt, &(pardisoData_->maxfct), 
		&(pardisoData_->mnum), &(pardisoData_->mtype), 
		&(pardisoData_->phase), &n,
		a, ia, ja, &(pardisoData_->idum), &(pardisoData_->nrhs), 
		pardisoData_->iparm, &(pardisoData_->msglvl), 
		b, x, &(pardisoData_->error));
	// Check for error
	if ( pardisoData_->error != 0 ) {
		printf ("\nERROR during solution: %d", pardisoData_->error);
	}
	// return
	return 0;
}

int destroyPardisoForSymmetricSystem(
		struct pardisoData* pardisoData_,
		int n, int* ia, int* ja) {
	//
	// Termination and release of memory
	// 
	// Release internal memory
	pardisoData_->phase = -1;
	PARDISO(pardisoData_->pt, &(pardisoData_->maxfct), &(pardisoData_->mnum), 
		&(pardisoData_->mtype), &(pardisoData_->phase),
		&n, &(pardisoData_->ddum), 
		ia, ja, &(pardisoData_->idum), &(pardisoData_->nrhs),
		pardisoData_->iparm, &(pardisoData_->msglvl), &(pardisoData_->ddum), 
		&(pardisoData_->ddum), &(pardisoData_->error));	
	// return
	return 0;
}

int allSolveSymmetricSystemWithPardiso(
		int n, int* ia, int* ja, double* a,
		double* b, double* x) {
	//
	// data structure
	//
	struct pardisoData pardisoData_;
	//
	// init
	//
	initPardisoInternalDataForSymmetricLinearSystem(&pardisoData_, n, ia, ja, a);
	//
	// solve
	//
	solveSymmetricSystemWithPardiso(&pardisoData_, n, ia, ja, a, b, x);
	//
	// destroy
	//
	destroyPardisoForSymmetricSystem(&pardisoData_, n, ia, ja);
	//
	// return
	//
	return 0;
}

void updateSecondStageInformationIfNeeded(
		struct DataSolverUpperSmoothing* user_data,
		Number* x_trial) {
	//
	// calc norm of x_trial - last_point_in_buffer
	//
	int n = user_data->totalDecisionVarsFirstStage;
	double normDiff = 0.0;
	for (int i=0; i<n; ++i) {
		double tmp = user_data->lastFirstStageTrial[i];
		tmp -= x_trial[i];
		tmp *= tmp;
		normDiff += tmp;
	}
	normDiff = sqrt(normDiff);
	//
	// verify if x_trial is a new point
	//
	if (normDiff > 0.000000001) {
		//
		// reset value and grad second stage buffer information
		//
		user_data->lastValuesPerScenario = 0.0;
		for (int  i=0; i<n; ++i) {
			user_data->lastGradSecondStage[i] = 0.0;
		}
		//
		// reset second order information second stage
		//
		if (user_data->useSecondOrderInfoSecondStage == true) {
			int size = user_data->totalPositionsHessSecondStage;
			for (int i=0; i<size; ++i) {
				user_data->valsHessSecondStage[i] = 0.0;
			}
		}
		// 
		// global indexes
		// 
		int globalIndexDecisionVars = 0;
		int globaIndexRhs = 0;
		int globalIndexMatrixT = 0;
		int globalIndexMatrixTt = 0;
		//
		// get references for the spaces of the right hand side and solution
		//
		double* b = user_data->rhsForNewtonPerScenario;
		double* sol = user_data->solForNewtonPerScenario;
		//
		// iterate scenarios
		//
		for (int s=0; s<(user_data->totalScenarios); ++s) {
			//
			// update flags inside main data structure
			//
			user_data->levelIter = 2;
			user_data->scenarioIter = s;
			//
			// get pointer to correct info for scenario
			//
			struct QpProbData* dataScenario = &(user_data->secondStageData[s]);
			int totDecisionVars = dataScenario->totalDecisionVars;
			int totalNzMatrixT = user_data->totalNonZerPerScenarioT[s];
			int innerCounterMatrixTtScen = 0;
			//
			// init right hand side vector with fixed rhs
			//
			for (int i=0; i<(dataScenario->totalEqualityConstraints); ++i) {
				dataScenario->rhsVec[i] = user_data->rhsPerScenario[globaIndexRhs + i];
			}
			//
			// add first stage information to rhs
			//
			for (int i=0; i<totalNzMatrixT; ++i) {
				int row = user_data->nzRowsPerScenarioT[globalIndexMatrixT + i];
				int col = user_data->nzColsPerScenarioT[globalIndexMatrixT + i];
				double val = user_data->nzValsPerScenarioT[globalIndexMatrixT + i];
				dataScenario->rhsVec[row-1] -= val * (x_trial[col-1]);
			}
			//
			// set pointers to the correct places
			//
			double optValScenIter = 0.0;
			double* optPrimalSolScenIter = &(user_data->optPrimalSolScen[globalIndexDecisionVars]);
			double* optMultEqConstraintsScenIter = &(user_data->optMultEqConstraintsScen[globaIndexRhs]);
			//
			// solve subproblem
			//
			solveMuTargetSubproblemWithIpopt(
				s, user_data,
				optPrimalSolScenIter,
				&optValScenIter,
				optMultEqConstraintsScenIter);
			//
			// to avoid diving by zero
			// 
			double MIN_VALUE = 0.00000001;
			shiftVectorToPositiveOrthantIfNeeded(
				totDecisionVars, optPrimalSolScenIter, MIN_VALUE);
			//
			// refine ipopt solution of mu-target subproblem with a damped newton method
			// NOTE: don't know how to configure ipopt tolerances to avoid this refinement
			// NOTE: but in any case, maybe it is good to have it separately
			//
			bool REFINE_IPOPT_SOL = true;
			if (REFINE_IPOPT_SOL == true) {
				refineMuTargetSolutionWithDampedNewtonMethod(
					s, user_data,
					optPrimalSolScenIter,
					optMultEqConstraintsScenIter);
				//
				// to avoid diving by zero
				// 
				shiftVectorToPositiveOrthantIfNeeded(
					totDecisionVars, optPrimalSolScenIter, MIN_VALUE);
			}
			//
			// increment total optimal value
			//
			user_data->lastValuesPerScenario += (user_data->probVecScenarios[s])*(optValScenIter);
			//
			// get pre-computed newton matrix in symmetric form
			// 
			int totalNewtonVars = user_data->sizesNewtonSystemsPerScenario[s];
			int* ia = user_data->rowIncrementsNewtonSystemsPerScenario[s];
			int* ja = user_data->colsNewtonSystemsPerScenario[s];
			double* a = user_data->valsNewtonSystemsPerScenario[s];
			int* positionsToUpdate = user_data->positionsOfDiagonalElementsToChangePerScenario[s];
			// 
			// update matrix to solve linear system of derivatives
			// 
			for (int i=0; i<(dataScenario->totalDecisionVars); ++i) {
				int index = positionsToUpdate[i];
				double tmp = (dataScenario->muTarget) / (optPrimalSolScenIter[i]);
				tmp = tmp / (optPrimalSolScenIter[i]);
				a[index] += tmp;
			}
			//
			// pardiso data
			//
			struct pardisoData pardisoData_;
			//
			// init pardiso and compute factorization
			// 
			initPardisoInternalDataForSymmetricLinearSystem(&pardisoData_, totalNewtonVars, ia, ja, a);
			//
			// part relative to primal vars is always zero
			//
			for (int j=0; j<totDecisionVars; ++j) {
				b[j] = 0.0;
			}
			//
			// iterate the coordinates of the first stage variable
			// 
			for (int i=0; i<n; ++i) {
				//
				// fill part relative to dual variables
				// 
				for (int j=totDecisionVars; j<totalNewtonVars; ++j) {
					b[j] = 0.0;
				}
				// 
				// fill right hand side vector
				// NOTE: theta_j = 0, varphi_j = 0
				// NOTE: phi_j = partial_j b (only non-zero part)
				// NOTE: partial_j b = (-) j-th column matrix T
				// NOTE: this is the reason we need T^t pre-computed
				// 
				while (innerCounterMatrixTtScen < totalNzMatrixT &&
						user_data->nzRowsPerScenarioTt[globalIndexMatrixTt] == i+1) {
					//
					// get data
					//
					int col = user_data->nzColsPerScenarioTt[globalIndexMatrixTt];
					double val = user_data->nzValsPerScenarioTt[globalIndexMatrixTt];
					//
					// set values
					//
					b[totDecisionVars+col-1] = - val;
					//
					// update counters
					//
					globalIndexMatrixTt += 1;
					innerCounterMatrixTtScen += 1;
				}
				//
				// calculate first order derivatives of primal-dual map
				// 
				solveSymmetricSystemWithPardiso(&pardisoData_, totalNewtonVars, ia, ja, a, b, sol);
				//
				// store solution for the computation of second order information
				// 
				if (user_data->useSecondOrderInfoSecondStage == true) {
					double* tmpBufferFirDer = user_data->derivativesPrimalPerScenarioWrtPcoord[i];
					for (int j=0; j<totDecisionVars; ++j) {
						tmpBufferFirDer[j] = sol[j];
					}
				}
				//
				// compute linear part of partial derivative of objective 
				//
				double partialDerivativeObjective = 0;
				for (int j=0; j<totDecisionVars; ++j) {
					double tmpVal = dataScenario->costVec[j];
					partialDerivativeObjective += (sol[j])*(tmpVal);
				}
				//
				// compute quadratic part of partial derivative of objective 
				//
				for (int j=0; j < (dataScenario->totalNzQuadCost) ; ++j) {
					int row = dataScenario->idxRowNzQuadCost[j];
					int col = dataScenario->idxColNzQuadCost[j];
					double val = dataScenario->valNzQuadCost[j];
					partialDerivativeObjective += 2*val*(sol[row-1])*(optPrimalSolScenIter[col-1]);
				}
				//
				// accumulate result for the gradient considering probability of scenario
				//
				double tmpVal = partialDerivativeObjective;
				tmpVal *= user_data->probVecScenarios[s];
				user_data->lastGradSecondStage[i] += tmpVal;
			}
			//
			// reset right hand side for computation of second order info
			// NOTE: Phi_kj is always zero and therefore the dual part is always zero
			//
			for (int j=totDecisionVars; j<totalNewtonVars; ++j) {
				b[j] = 0.0;
			}
			//
			// compute hessian of optimal value of regularized recourse problem
			// 
			if (user_data->useSecondOrderInfoSecondStage == true) {
				int globalIndexSecondOrderInfo = 0;
				double localContrib = 0;
				for (int row=1; row<=n; ++row) {
					for (int col=row; col<=n; ++col) {
						//
						// get pointers to stored first order derivatives
						//
						double* firtDer1 = user_data->derivativesPrimalPerScenarioWrtPcoord[row-1];
						double* firtDer2 = user_data->derivativesPrimalPerScenarioWrtPcoord[col-1];
						// 
						// compute right hand side for second order info
						// NOTE: only non-zero, but diagonal, term is partial_p_k M
						// 
						for (int j=0; j<totDecisionVars; ++j) {
							// 
							// compute diagonal of \partial_p_k M
							// 
							double tmp1 = firtDer1[j];
							double tmp2 = firtDer2[j];
							double diagTerm = 2*(dataScenario->muTarget);
							diagTerm = diagTerm / optPrimalSolScenIter[j];
							diagTerm = diagTerm / optPrimalSolScenIter[j];
							diagTerm = diagTerm / optPrimalSolScenIter[j];
							diagTerm = diagTerm * tmp1 * tmp2;
							//
							// set rhs
							//	
							b[j] = diagTerm;
						}
						//
						// compute second order derivative of solution mappings
						//
						solveSymmetricSystemWithPardiso(&pardisoData_, 
							totalNewtonVars, ia, ja, a, b, sol);
						//
						// init local contribution
						//
						localContrib = 0.0;
						//
						// calc local contrib relative to the linear part of the objective
						//
						for (int j=0; j<totDecisionVars; ++j) {
							localContrib += (dataScenario->costVec[j])*(sol[j]);
						}
						// 
						// calc local contrib relative to the quadratic part of the objective
						//
						for (int j=0; j < (dataScenario->totalNzQuadCost) ; ++j) {
							int rowNzQuad = dataScenario->idxRowNzQuadCost[j];
							int colNzQuad = dataScenario->idxColNzQuadCost[j];
							double valNzQuad = dataScenario->valNzQuadCost[j];
							localContrib += valNzQuad * (sol[rowNzQuad-1]) * (optPrimalSolScenIter[colNzQuad-1]);
							localContrib += valNzQuad * (sol[colNzQuad-1]) * (optPrimalSolScenIter[rowNzQuad-1]);
							localContrib += valNzQuad * firtDer1[rowNzQuad-1] * firtDer2[colNzQuad-1];
							localContrib += valNzQuad * firtDer1[colNzQuad-1] * firtDer2[rowNzQuad-1];
						}
						// 
						// set second order info second stage
						// NOTE: here we just append the second order information
						// NOTE: because it was initialized with zero in the beginning
						//
						localContrib *= (user_data->probVecScenarios[s]);
						user_data->valsHessSecondStage[globalIndexSecondOrderInfo] += localContrib;
						//
						// update global index
						//
						globalIndexSecondOrderInfo += 1;
					}
				}
			}
			//
			// release memory inside pardiso
			//
			destroyPardisoForSymmetricSystem(&pardisoData_, totalNewtonVars, ia, ja);
			//
			// restore newton matrix to initial state
			//
			for (int i=0; i<(dataScenario->totalDecisionVars); ++i) {
				int index = positionsToUpdate[i];
				double tmp = (dataScenario->muTarget) / (optPrimalSolScenIter[i]);
				tmp = tmp / (optPrimalSolScenIter[i]);
				a[index] -= tmp;
			}
			//
			// update global indexes
			//
			globaIndexRhs += dataScenario->totalEqualityConstraints;
			globalIndexMatrixT += user_data->totalNonZerPerScenarioT[s];
			globalIndexDecisionVars += dataScenario->totalDecisionVars;
		}
		// 
		// update point in buffer
		// 
		for (int i=0; i<n; ++i) {
			user_data->lastFirstStageTrial[i] = x_trial[i];
		}
		//
		// move flags inside main data structure
		// back to first stage configuration
		// 
		user_data->levelIter = 1;
		user_data->scenarioIter = -1;
	}
}

// Function Implementations 
Bool eval_f_quad_prob(Index n, Number* x, Bool new_x, Number* obj_value, UserDataPtr user_data_all) {
	// cast
	struct DataSolverUpperSmoothing* user_data_all_cast = (DataSolverUpperSmoothing*) user_data_all;
	// get specific data
	struct QpProbData* data = NULL;
	if (user_data_all_cast->levelIter == 1) {
		data = &(user_data_all_cast->firstStageData);
		updateSecondStageInformationIfNeeded(
			user_data_all_cast, x);
	} else {
		int scenario = user_data_all_cast->scenarioIter;
		data = &(user_data_all_cast->secondStageData[scenario]);
	}
	// init	
	Number objValTmp = 0.0;
	// calc linear part
	for (int i=0; i < n; ++i) {
		double tmpVal = (data->costVec[i])*(x[i]);
		objValTmp += tmpVal;
	}

	// calc quadratic part
	int totNzQuad = data->totalNzQuadCost;
	for (int i=0; i<totNzQuad; ++i) {
		int row = data->idxRowNzQuadCost[i];
		int col = data->idxColNzQuadCost[i];
		double val = data->valNzQuadCost[i];
		Number tmpVal = val*(x[row-1])*(x[col-1]);
		objValTmp += tmpVal;
	}
	// correction only at the first level
	if (user_data_all_cast->levelIter == 1) {
		// add gradient due to recourse function of scenarios
		objValTmp += user_data_all_cast->lastValuesPerScenario;
	}
	// set
	*obj_value = objValTmp;
	// return
	return TRUE;
}

Bool eval_grad_f_quad_prob(Index n, Number* x, Bool new_x, Number* grad_f, UserDataPtr user_data_all) {
	// cast
	struct DataSolverUpperSmoothing* user_data_all_cast = (DataSolverUpperSmoothing*) user_data_all;
	// get specific data
	struct QpProbData* data = NULL;
	if (user_data_all_cast->levelIter == 1) {
		data = &(user_data_all_cast->firstStageData);
		updateSecondStageInformationIfNeeded(
			user_data_all_cast, x);
	} else {
		int scenario = user_data_all_cast->scenarioIter;
		data = &(user_data_all_cast->secondStageData[scenario]);
	}
	// init with linear part
	for (int i=0; i < n; ++i) {
		grad_f[i] = data->costVec[i];
	}
	// add quadratic contribution
	int totNzQuad = data->totalNzQuadCost;
	for (int i=0; i<totNzQuad; ++i) {
		int row = data->idxRowNzQuadCost[i];
		int col = data->idxColNzQuadCost[i];
		double val = data->valNzQuadCost[i];
		grad_f[row-1] += 2*val*x[col-1];
	}
	// correction only at the first level
	if (user_data_all_cast->levelIter == 1) {
		// add gradient due to recourse function of scenarios
		for (int i=0; i < n; ++i) {
			grad_f[i] += user_data_all_cast->lastGradSecondStage[i];
		}
	}
	// return
	return TRUE;
}

Bool eval_g_quad_prob(Index n, Number* x, Bool new_x, Index m, Number* g, UserDataPtr user_data_all) {
	// cast
	struct DataSolverUpperSmoothing* user_data_all_cast = (DataSolverUpperSmoothing*) user_data_all;
	// get specific data
	struct QpProbData* data = NULL;
	if (user_data_all_cast->levelIter == 1) {
		data = &(user_data_all_cast->firstStageData);
	} else {
		int scenario = user_data_all_cast->scenarioIter;
		data = &(user_data_all_cast->secondStageData[scenario]);
	}
	// init
	for (int i=0; i < m; ++i) {
		g[i] = - data->rhsVec[i];
	}
	// evaluate constraints
	int totNzMatrixA = data->totalNzMatrixA;
	for (int i=0; i < totNzMatrixA; ++i) {
		int row = data->idxRowNzMatrixA[i];
		int col = data->idxColNzMatrixA[i];
		double val = data->valNzMatrixA[i];
		g[row-1] += val * (x[col-1]);
	}
	// return
	return TRUE;
}

Bool eval_jac_g_quad_prob(Index n, Number *x, Bool new_x, Index m, Index nele_jac, 
		Index *iRow, Index *jCol, Number *values, UserDataPtr user_data_all) {
	// cast
	struct DataSolverUpperSmoothing* user_data_all_cast = (DataSolverUpperSmoothing*) user_data_all;
	// get specific data
	struct QpProbData* data = NULL;
	if (user_data_all_cast->levelIter == 1) {
		data = &(user_data_all_cast->firstStageData);
	} else {
		int scenario = user_data_all_cast->scenarioIter;
		data = &(user_data_all_cast->secondStageData[scenario]);
	}
	// sparsity pattern
	if (values == NULL) {
		int totNzMatrixA = data->totalNzMatrixA;
		for (int i=0; i < totNzMatrixA; ++i) {
			iRow[i] = data->idxRowNzMatrixA[i] - 1;
			jCol[i] = data->idxColNzMatrixA[i] - 1;
		}
	// values
	} else {
		int totNzMatrixA = data->totalNzMatrixA;
		for (int i=0; i < totNzMatrixA; ++i) { 
			values[i] = data->valNzMatrixA[i];
		}
	}
	// return
	return TRUE;
}

Bool eval_h_quad_prob(Index n, Number *x, Bool new_x, Number obj_factor, Index m, 
		Number *lambda, Bool new_lambda, Index nele_hess, 
		Index *iRow, Index *jCol, Number *values, UserDataPtr user_data_all) {
	// cast
	struct DataSolverUpperSmoothing* user_data_all_cast = (DataSolverUpperSmoothing*) user_data_all;
	// get specific data
	struct QpProbData* data = NULL;
	if (user_data_all_cast->levelIter == 1 &&
			user_data_all_cast->useSecondOrderInfoSecondStage == true) {
		//
		// here we have the quadratic first stage quadratic perturbation
		// NOTE: but this is currently not being used
		// NOTE: so we just fill the second stage terms
		//
		// TODO : data = &(user_data_all_cast->firstStageData);
		int size = user_data_all_cast->totalPositionsHessSecondStage;
		// sparsity pattern
		if (values == NULL) {
			for (int i=0; i<size; ++i) {
				iRow[i] = user_data_all_cast->idxRowHessSecondStage[i] - 1;
				jCol[i] = user_data_all_cast->idxColHessSecondStage[i] - 1;
			}
		// fill values
		} else {
			for (int i=0; i<size; ++i) {
				double val = user_data_all_cast->valsHessSecondStage[i];
				values[i] = obj_factor * val;
			}
		}
	} else {
		int scenario = user_data_all_cast->scenarioIter;
		data = &(user_data_all_cast->secondStageData[scenario]);
		// sparsity pattern
		if (values == NULL) {
			for (int i=0; i < data->totalNzSymQuadCost; ++i) {
				iRow[i] = data->idxRowNzSymQuadCost[i] - 1;
				jCol[i] = data->idxColNzSymQuadCost[i] - 1;
			}
		// fill values
		} else {
			for (int i=0; i < data->totalNzSymQuadCost; ++i) {
				double tmpVal = data->valNzSymQuadCost[i];
				values[i] = obj_factor * 2 * tmpVal;
			}
		}
	}
	// return
	return TRUE;
}

int getStartingPointForIpoptWithCplex(
		struct QpProbData* data,
		double* initialPrimalIterate,
		double* initialDualIterate) {
	// 
	// get dimensions of quadratic problem
	// 
	int totDecVarIter = data->totalDecisionVars;
	int totEqConstrIter = data->totalEqualityConstraints;
	int totalNzMatrixA = data->totalNzMatrixA;
	//
	// try, because this is life
	//
	try {
		//
		// build model
		//
		IloEnv env;
		IloModel model(env);
		//
		// create var
		//
		IloNumVarArray x(env);
		for(int i=0; i<totDecVarIter; ++i) {
			x.add(IloNumVar(env));
		}
		//
		// create objective for cplex solver
		//
		IloExpr objExpr(env);
		for (int i=0; i<totDecVarIter; ++i) {
			objExpr += (data->costVec[i])*(x[i]);
		}
		for (int i=0; i<(data->totalNzQuadCost); ++i) {
			int row = data->idxRowNzQuadCost[i];
			int col = data->idxColNzQuadCost[i];
			double val = data->valNzQuadCost[i];
			objExpr += val*(x[row-1])*(x[col-1]);
		}
		//
		// set cplex objective
		//
		IloObjective obj = IloMinimize(env, objExpr);
		model.add(obj);
		//
		// set non-negativity constraints for variables of second stage problems
		//
		for (int i=0; i<totDecVarIter; ++i) {
			model.add(x[i] >= 0);
		}
		//
		// create equality constraints
		//
		int globalIndexScenarioMatrixA = 0;
		IloRangeArray constraints(env);
		for(int i=0; i<totEqConstrIter; ++i) {
			//
			// compute expression
			//
			IloExpr linExpr(env);
			int counter = 0;
			while(counter < totalNzMatrixA && data->idxRowNzMatrixA[globalIndexScenarioMatrixA] == i+1) {
				int rowIter = data->idxRowNzMatrixA[globalIndexScenarioMatrixA];
				int colIter = data->idxColNzMatrixA[globalIndexScenarioMatrixA];
				double valIter = data->valNzMatrixA[globalIndexScenarioMatrixA];
				linExpr += valIter*x[colIter-1];
				globalIndexScenarioMatrixA++;
				counter++;
			}
			//
			// add constraint
			//
			constraints.add(linExpr == data->rhsVec[i]);
			//
			//
			//
			linExpr.end();
		}
		//
		// add constraints
		//
		model.add(constraints);
		//
		// solve scenario lp
		//
		IloCplex cplex(model);
		cplex.setParam(IloCplex::RootAlg, IloCplex::Param::Barrier::Algorithm);
		cplex.setParam(IloCplex::Param::Threads, 1);
		cplex.setParam(IloCplex::Param::Barrier::ConvergeTol, 0.0000000001);
		IloBool isFeasible = cplex.solve();
		//
		// get primal solution
		//
		for(int i=0; i<totDecVarIter; ++i) {
			initialPrimalIterate[i] = cplex.getValue(x[i]);
		}
		//
		// get dual solution
		// 
		for(int i=0; i<totEqConstrIter; ++i) {
			initialDualIterate[i] = cplex.getValue(constraints[i]);
		}
		// end
		objExpr.end();
		obj.end();
		x.end();
		constraints.end();
		model.end();
		cplex.end();
		env.end();
	} catch (...) {
		//
		// get primal solution
		//
		for(int i=0; i<totDecVarIter; ++i) {
			initialPrimalIterate[i] = 0.0;
		}
		//
		// get dual solution
		// 
		for(int i=0; i<totEqConstrIter; ++i) {
			initialDualIterate[i] = 0.0;
		}
	}
}

void computeGradientOfLagrangian(
		double* gradientOfLagrangian,
		struct QpProbData* data,
		double* primalSolIter,
		double* multiplierIter) {
	//
	// init primal part with linear cost
	//
	int totalDecisionVars = data->totalDecisionVars;
	int totalEqualityConstraints = data->totalEqualityConstraints;
	int totNewtonVars = totalDecisionVars + totalEqualityConstraints;
	for (int i=0; i<totalDecisionVars; ++i) {
		gradientOfLagrangian[i] = data->costVec[i];
	}
	//
	// init dual part with -b
	//
	for (int i=totalDecisionVars;i<totNewtonVars; ++i) {
		gradientOfLagrangian[i] = (-1.0) * (data->rhsVec[i - totalDecisionVars]);
	}
	//
	// fill part due to quadratic cost
	//
	int totalNzQuadCost = data->totalNzQuadCost;
	for (int i=0; i<totalNzQuadCost; ++i) {
		int row = data->idxRowNzQuadCost[i];
		int col = data->idxColNzQuadCost[i];
		double val = data->valNzQuadCost[i];
		gradientOfLagrangian[row-1] += 2 * val * primalSolIter[col-1];
	}
	//
	// fill part due to x >= 0
	//
	for (int i=0; i<totalDecisionVars; ++i) {
		gradientOfLagrangian[i] -= (data->muTarget) / primalSolIter[i];
	}
	//
	// fill part due to +A^t lambda
	//
	int totalNzMatrixA = data->totalNzMatrixA;
	for (int i=0;i<totalNzMatrixA; ++i) {
		int row = data->idxRowNzMatrixA[i];
		int col = data->idxColNzMatrixA[i];
		double val = data->valNzMatrixA[i];
		gradientOfLagrangian[col-1] += val * (multiplierIter[row-1]);
	}
	//
	// fill part due to : Ax
	//
	for (int i=0;i<totalNzMatrixA; ++i) {
		int row = data->idxRowNzMatrixA[i];
		int col = data->idxColNzMatrixA[i];
		double val = data->valNzMatrixA[i];
		gradientOfLagrangian[totalDecisionVars+row-1] += val * (primalSolIter[col-1]);
	}
}

int refineMuTargetSolutionWithDampedNewtonMethod(
		int scenario,
		struct DataSolverUpperSmoothing* user_data,
		double* initialPrimalIterate,
		double* initialDualIterate) {
	//
	// get scenario data
	//
	struct QpProbData* data = &(user_data->secondStageData[scenario]);
	//
	// constants
	//
	double KKT_ERROR_TOLERANCE = 0.00000000001;
	double ALLOWED_PROP_DECREASE_PRIMAL_ITERATE = 0.01;
	double DECREASE_IN_DAMPING_BACKTRACKING = 0.9;
	double MAX_ITER_STEP_SIZE_BACKTRACK = 200;
	double STEP_SIZE_REDUCTION = 0.95;
	bool PRINT_KKT_ERROR = false;
	int MAX_ITER = data->totalNewtonRefinementStepsAfterIpopt;
	//
	// init flags
	//
	int totalDecisionVars = data->totalDecisionVars;
	int totalEqualityConstraints = data->totalEqualityConstraints;
	int totNewtonVars = totalDecisionVars + totalEqualityConstraints;
	double kktError = 1000000.0;
	int totalIter = 0;
	//
	// allocate vectors
	//
	double* primalSolIter = user_data->primalSolIter;
	double* multiplierIter = user_data->multiplierIter; 
	double* gradientOfLagrangianIter = user_data->gradientOfLagrangianIter;
	double* gradLagrTrialPoint = user_data->gradLagrTrialPoint;
	double* newtonDirectionIter = user_data->newtonDirectionIter;
	double* primalSolIterDamp = user_data->primalSolIterDamp;
	double* multIterDamp = user_data->multIterDamp;
	//
	// init
	//
	for (int i=0; i<totalDecisionVars; ++i) {
		primalSolIter[i] = initialPrimalIterate[i];
	}
	for (int i=0; i<totalEqualityConstraints; ++i) {
		multiplierIter[i] = initialDualIterate[i];
	}
	// 
	// alloc data to send to pardiso mkl
	//
	int totNzNewtonMatrix = user_data->sizesNewtonSystemsPerScenario[scenario];
	MKL_INT* rowIncrementNewtonMatrix = user_data->rowIncrementsNewtonSystemsPerScenario[scenario];
	MKL_INT* colIdxNewtonMatrix = user_data->colsNewtonSystemsPerScenario[scenario];
	double* valsNewtonMatrix = user_data->valsNewtonSystemsPerScenario[scenario];
	int* positionsToChange = user_data->positionsOfDiagonalElementsToChangePerScenario[scenario];
	// 
	// iterate
	//
	while (kktError > KKT_ERROR_TOLERANCE && totalIter < MAX_ITER) {
		//
		// compute gradient of lagrangian
		//
		computeGradientOfLagrangian(
			gradientOfLagrangianIter,
			data,
			primalSolIter,
			multiplierIter);
		//
		// change sign
		//
		for (int i=0; i<totNewtonVars; ++i) {
			gradientOfLagrangianIter[i] *= -1.0;
		}
		//
		// update KKT error
		//
		kktError = normVec(gradientOfLagrangianIter, totNewtonVars);
		//
		// add to newton matrix the X^(-2.0) contribution
		//
		for (int i=0;i<totalDecisionVars; ++i) {
			int index = positionsToChange[i];
			valsNewtonMatrix[index] += (data->muTarget) / ((primalSolIter[i])*(primalSolIter[i]));
		}
		//
		// compute newton direction with pardiso mkl
		//
		allSolveSymmetricSystemWithPardiso(
			totNewtonVars, rowIncrementNewtonMatrix, 
			colIdxNewtonMatrix, valsNewtonMatrix, 
			gradientOfLagrangianIter, newtonDirectionIter);
		//
		// remove from newton matrix the X^(-2.0) contribution
		//
		for (int i=0;i<totalDecisionVars; ++i) {
			int index = positionsToChange[i];
			valsNewtonMatrix[index] -= (data->muTarget) / ((primalSolIter[i])*(primalSolIter[i]));
		}
		//
		// get step size that does not violate proportional decrease
		//
		double stepSize = 1.0;
		for (int i=0; i<totalDecisionVars; ++i) {
			if (newtonDirectionIter[i] < 0) {
				double tmp = ALLOWED_PROP_DECREASE_PRIMAL_ITERATE - 1;
				tmp = tmp / newtonDirectionIter[i];
				tmp = tmp * primalSolIter[i];
				if (tmp < stepSize) {
					stepSize = tmp;
				}
			}
		}
		//
		// linesearch, backtracking part of damped newton method
		//
		bool foundStepSize = false;
		int totalIterBacktrack = 0;
		while (foundStepSize == false && totalIterBacktrack <= MAX_ITER_STEP_SIZE_BACKTRACK) {
			//
			// compute primal trial point
			//
			for (int i=0; i<totalDecisionVars; ++i) {
				primalSolIterDamp[i] = 	primalSolIter[i] + stepSize*(newtonDirectionIter[i]);
			}
			//
			// compute dual trial point
			//
			for (int i=0; i<totalEqualityConstraints; ++i) {
				multIterDamp[i] = multiplierIter[i] + stepSize*(newtonDirectionIter[totalDecisionVars+i]);
			}
			//
			// compute gradient of lagrangian at trial point
			//
			computeGradientOfLagrangian(
				gradientOfLagrangianIter,
				data,
				primalSolIter,
				multiplierIter);
			// 
			// test decrease newton system
			// 
			double normTrial = normVec(gradLagrTrialPoint, totNewtonVars);
			if (normTrial <= DECREASE_IN_DAMPING_BACKTRACKING*kktError) {
				foundStepSize = true;
			} else {
				stepSize *= STEP_SIZE_REDUCTION;
			}
			//
			// update
			//
			totalIterBacktrack += 1;
		}
		// 
		// kkt error output
		//
		if (PRINT_KKT_ERROR == true) {
			printf(" \n kktError : %f \n", kktError);
			printf(" stepSize : %f \n", stepSize);
		}
		//
		// update primal iterate
		// 
		for (int i=0; i<totalDecisionVars; ++i) {
			primalSolIter[i] += stepSize * newtonDirectionIter[i];
		}	
		//
		// update dual iterate
		//
		for (int i=0; i<totalEqualityConstraints; ++i) {
			multiplierIter[i] += stepSize * newtonDirectionIter[i + totalDecisionVars];
		}
		//
		// update counter
		//
		totalIter += 1;
	}
	//
	// copy primal sol
	//
	for (int i=0; i<totalDecisionVars; ++i) {
		initialPrimalIterate[i] = primalSolIter[i];
	}
	//
	// copy dual sol with inverse signal
	//
	for (int i=0; i<totalEqualityConstraints; ++i) {
		initialDualIterate[i] = multiplierIter[i];
	}
	//
	// return
	//
	return 0;
}

int solveMuTargetSubproblemWithIpopt(
		int scenario,
		struct DataSolverUpperSmoothing* user_data,
		double* initialPrimalIterate,
		double* optVal,
		double* initialDualIterate) {
	// get specific data
	struct QpProbData* user_data_iter = NULL;
	if (scenario == -1) {
		user_data_iter = &(user_data->firstStageData);
	} else {
		user_data_iter = &(user_data->secondStageData[scenario]);
	}
	//
	// get good and feasible starting point for ipopt
	//
	getStartingPointForIpoptWithCplex(
			user_data_iter,
			initialPrimalIterate,
			initialDualIterate);
	// number of variables 
	Index n=-1;
	 // number of constraints 
	Index m=-1;
	// lower bounds on x               
	Number* x_L = NULL;
	// upper bounds on x    
	Number* x_U = NULL;
	// lower bounds on g 
	Number* g_L = NULL;
	// upper bounds on g 
	Number* g_U = NULL;
	// IpoptProblem 
	IpoptProblem nlp = NULL;
	// Solve return code 
	enum ApplicationReturnStatus status;
	// lower bound multipliers at the solution 
	Number* mult_x_L = NULL;
	// upper bound multipliers at the solution
	Number* mult_x_U = NULL;
	// objective value  
	Number obj = 0.0;
	// generic counter
	Index i;	
	// Number of decision variables
	n = user_data_iter->totalDecisionVars;
	// Number of constraints
	m = user_data_iter->totalEqualityConstraints;
	// Number of nonzeros in the Jacobian of the constraints 
	Index nele_jac = user_data_iter->totalNzMatrixA;
	// Number of nonzeros in the Hessian of the Lagrangian (lower or upper triangual part only) 
	Index nele_hess = 0;
	if (scenario == -1 && user_data->useSecondOrderInfoSecondStage == true) {
		// 
		// in this case we consider a full dense symmetric hessian 
		// 
		nele_hess = n*(n+1) / 2;
	} else if (scenario != -1 && user_data_iter->useSecondOrderInfo == true) {
		// 
		// for scenario problems we just put the sparse quadratic term
		//
		nele_hess = user_data_iter->totalNzSymQuadCost;
	}

	// indexing style for matrices 
	Index index_style = 0;
	// set the number of variables and allocate space for the bounds 
	x_L = (Number*)malloc(sizeof(Number)*n);
	x_U = (Number*)malloc(sizeof(Number)*n);
	// set the values for the variable bounds 
	for (i=0; i<n; i++) {
		x_L[i] = 0.0;
		x_U[i] = 2e19;
	}
	// set the number of constraints and allocate space for the bounds 
	if (m > 0) {
		g_L = (Number*)malloc(sizeof(Number)*m);
		g_U = (Number*)malloc(sizeof(Number)*m);
		// set the values of the constraint bounds 
		for (int i=0; i<m; ++i) {
			g_L[i] = 0.0;
			g_U[i] = 0.0;
		}
	}
	// create the IpoptProblem
	nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
		index_style, &eval_f_quad_prob, &eval_g_quad_prob, &eval_grad_f_quad_prob,
		&eval_jac_g_quad_prob, &eval_h_quad_prob);
	// We can free the memory now - the values for the bounds have been
	// copied internally in CreateIpoptProblem
	free(x_L);
	free(x_U);
	if (m > 0) {
		free(g_L);
		free(g_U);
	}
	// Set some options.  Note the following ones are only examples,
	// they might not be suitable for your problem.
	AddIpoptIntOption(nlp, "print_level", 0);
	AddIpoptIntOption(nlp, "max_iter", user_data_iter->maxIterIpopt);
	AddIpoptNumOption(nlp, "mu_target", user_data_iter->muTarget);
	AddIpoptStrOption(nlp, "warm_start_init_point", "yes");
	// The following option reduce the automatic modification of the
	// starting point done my Ipopt. 
	AddIpoptNumOption(nlp, "bound_push", 1e-5);
	AddIpoptNumOption(nlp, "bound_frac", 1e-5);
	// for the master problem no hessian
	if (user_data_iter->useSecondOrderInfo == false) {
		AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");		
	} else {
		AddIpoptStrOption(nlp, "hessian_approximation", "exact");		
	}
	// allocate space to store the bound multipliers at the solution
	mult_x_L = (Number*)malloc(sizeof(Number)*n);
	mult_x_U = (Number*)malloc(sizeof(Number)*n);
	// solve the problem
	status = IpoptSolve(nlp, initialPrimalIterate, NULL, &obj, initialDualIterate, mult_x_L, mult_x_U, user_data);
	// copy optimal value
	*optVal = obj;
	// free 
	FreeIpoptProblem(nlp);
	free(mult_x_L);
	free(mult_x_U);
	// return
	return (int)status;
}

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
		double* nzValsQuadPerturbationPerScenario_) {
	//
	// the gross data
	//
	data_->useSecondOrderInfoSecondStage = useSecondOrderInfoMasterProblem;
	data_->epsilon = epsilon_;
	data_->tikhonovMu = tikhonovMu_;
	data_->totalDecisionVarsFirstStage = totalDecisionVarsFirstStage_;
	data_->totalScenarios = totalScenarios_;
	data_->totalDecisionVarsPerScenario = totalDecisionVarsPerScenario_;
	data_->totalEqConsFirstStage = totalEqConsFirstStage_;
	data_->rhsFirstStage = rhsFirstStage_;
	data_->totalEqConstrPerScenario = totalEqConstrPerScenario_;
	data_->costFirstStage = costFirstStage_;
	data_->probVecScenarios = probVecScenarios_;
	data_->costPerScenario = costPerScenario_;
	data_->rhsPerScenario = rhsPerScenario_;
	data_->totalNonZeroFirstStageW = totalNonZeroFirstStageW_;
	data_->nzRowsFirstStageW = nzRowsFirstStageW_;
	data_->nzColsFirstStageW = nzColsFirstStageW_;
	data_->nzValsFirstStageW = nzValsFirstStageW_;
	data_->totalNonZeroPerScenarioW = totalNonZeroPerScenarioW_;
	data_->nzRowsPerScenarioW = nzRowsPerScenarioW_;
	data_->nzColsPerScenarioW = nzColsPerScenarioW_;
	data_->nzValsPerScenarioW = nzValsPerScenarioW_;
	data_->nzRowsPerScenarioWt = nzRowsPerScenarioWt_;
	data_->nzColsPerScenarioWt = nzColsPerScenarioWt_;
	data_->nzValsPerScenarioWt = nzValsPerScenarioWt_;
	data_->totalNonZerPerScenarioT = totalNonZerPerScenarioT_;
	data_->nzRowsPerScenarioT = nzRowsPerScenarioT_;
	data_->nzColsPerScenarioT = nzColsPerScenarioT_;
	data_->nzValsPerScenarioT = nzValsPerScenarioT_;
	data_->nzRowsPerScenarioTt = nzRowsPerScenarioTt_;
	data_->nzColsPerScenarioTt = nzColsPerScenarioTt_;
	data_->nzValsPerScenarioTt = nzValsPerScenarioTt_;
	data_->totalNzQuadPerturbation = totalNzQuadPerturbation_;
	data_->nzRowsQuadPerturbation = nzRowsQuadPerturbation_;
	data_->nzColsQuadPerturbation = nzColsQuadPerturbation_;
	data_->nzValsQuadPerturbation = nzValsQuadPerturbation_;
	data_->totalNzQuadPerturbationPerScenario = totalNzQuadPerturbationPerScenario_;
	data_->nzRowsQuadPerturbationPerScenario = nzRowsQuadPerturbationPerScenario_;
	data_->nzColsQuadPerturbationPerScenario = nzColsQuadPerturbationPerScenario_;
	data_->nzValsQuadPerturbationPerScenario = nzValsQuadPerturbationPerScenario_;
	//
	// organize first stage data
	// 
	data_->firstStageData.useSecondOrderInfo = useSecondOrderInfoMasterProblem;
	data_->firstStageData.totalNewtonRefinementStepsAfterIpopt = 0;
	data_->firstStageData.maxIterIpopt = maxIterPerSubProblem;
	data_->firstStageData.muTarget = 0;
	data_->firstStageData.tikhonovMu = 0;
	data_->firstStageData.totalDecisionVars = totalDecisionVarsFirstStage_;
	data_->firstStageData.totalEqualityConstraints = totalEqConsFirstStage_;
	data_->firstStageData.costVec = costFirstStage_;
	data_->firstStageData.totalNzQuadCost = totalNzQuadPerturbation_;
	data_->firstStageData.idxRowNzQuadCost = nzRowsQuadPerturbation_;
	data_->firstStageData.idxColNzQuadCost = nzColsQuadPerturbation_;
	data_->firstStageData.valNzQuadCost = nzValsQuadPerturbation_;
	data_->firstStageData.totalNzSymQuadCost = 0;
	data_->firstStageData.idxRowNzSymQuadCost = NULL;
	data_->firstStageData.idxColNzSymQuadCost = NULL;
	data_->firstStageData.valNzSymQuadCost = NULL;
	data_->firstStageData.totalNzMatrixA = totalNonZeroFirstStageW_;
	data_->firstStageData.idxRowNzMatrixA = nzRowsFirstStageW_;
	data_->firstStageData.idxColNzMatrixA = nzColsFirstStageW_;
	data_->firstStageData.valNzMatrixA = nzValsFirstStageW_;
	data_->firstStageData.idxRowNzMatrixAt = NULL;
	data_->firstStageData.idxColNzMatrixAt = NULL;
	data_->firstStageData.valNzMatrixAt = NULL;
	data_->firstStageData.rhsVec = rhsFirstStage_;
	//
	// alloc data for symmetric matrix
	// 
	data_->firstStageData.idxRowNzSymQuadCost = (int*) malloc( sizeof(int) * totalNzQuadPerturbation_);
	data_->firstStageData.idxColNzSymQuadCost = (int*) malloc( sizeof(int) * totalNzQuadPerturbation_);
	data_->firstStageData.valNzSymQuadCost = (double*) malloc( sizeof(double) * totalNzQuadPerturbation_);
	//
	// init symmetric data
	// 
	for (int i=0; i < totalNzQuadPerturbation_; ++i) {
		//
		// get data
		//
		int row = nzRowsQuadPerturbation_[i];
		int col = nzColsQuadPerturbation_[i];
		double val = nzValsQuadPerturbation_[i];
		//
		// test if is upper triangle
		//
		if (col >= row) {
			data_->firstStageData.idxRowNzSymQuadCost[data_->firstStageData.totalNzSymQuadCost] = row;
			data_->firstStageData.idxColNzSymQuadCost[data_->firstStageData.totalNzSymQuadCost] = col;
			data_->firstStageData.valNzSymQuadCost[data_->firstStageData.totalNzSymQuadCost] = val;
			data_->firstStageData.totalNzSymQuadCost += 1;
		}
	}
	//
	// init global indexes
	//
	int globalIndexCostScenario = 0;
	int globalIndexQuadPerturbation = 0;
	int globalIndexMatrixA = 0;
	//
	// global counters
	//
	int totalDecisionVarsSecondStage = 0;
	int totalEqConstraintsSecondStage = 0;
	// 
	// iterate scenarios
	// 
	for (int s=0; s<totalScenarios_; ++s) {	
		//
		// set data
		//	
		struct QpProbData tmpDataQpScenario;
		tmpDataQpScenario.useSecondOrderInfo = true;
		tmpDataQpScenario.maxIterIpopt = maxIterPerSubProblem;
		tmpDataQpScenario.totalNewtonRefinementStepsAfterIpopt = 
			totalNewtonRefinementStepsAfterIpopt_;
		tmpDataQpScenario.muTarget = epsilon_;
		tmpDataQpScenario.tikhonovMu = tikhonovMu_;
		tmpDataQpScenario.totalDecisionVars = totalDecisionVarsPerScenario_[s];
		tmpDataQpScenario.totalEqualityConstraints = totalEqConstrPerScenario_[s];
		tmpDataQpScenario.costVec = &(costPerScenario_[globalIndexCostScenario]);
		tmpDataQpScenario.totalNzQuadCost = totalNzQuadPerturbationPerScenario_[s];
		tmpDataQpScenario.idxRowNzQuadCost = &(nzRowsQuadPerturbationPerScenario_[globalIndexQuadPerturbation]);
		tmpDataQpScenario.idxColNzQuadCost = &(nzColsQuadPerturbationPerScenario_[globalIndexQuadPerturbation]);
		tmpDataQpScenario.valNzQuadCost = &(nzValsQuadPerturbationPerScenario_[globalIndexQuadPerturbation]);
		tmpDataQpScenario.totalNzMatrixA = totalNonZeroPerScenarioW_[s];
		tmpDataQpScenario.idxRowNzMatrixA = &(nzRowsPerScenarioW_[globalIndexMatrixA]);
		tmpDataQpScenario.idxColNzMatrixA = &(nzColsPerScenarioW_[globalIndexMatrixA]);
		tmpDataQpScenario.valNzMatrixA = &(nzValsPerScenarioW_[globalIndexMatrixA]);
		tmpDataQpScenario.idxRowNzMatrixAt = &(nzRowsPerScenarioWt_[globalIndexMatrixA]);
		tmpDataQpScenario.idxColNzMatrixAt = &(nzColsPerScenarioWt_[globalIndexMatrixA]);
		tmpDataQpScenario.valNzMatrixAt = &(nzValsPerScenarioWt_[globalIndexMatrixA]);
		// 
		// update counters
		// 
		totalDecisionVarsSecondStage += tmpDataQpScenario.totalDecisionVars;
		totalEqConstraintsSecondStage += tmpDataQpScenario.totalEqualityConstraints;
		//
		// rhs set at the iteration
		// 
		int rhsSize = tmpDataQpScenario.totalEqualityConstraints;
		tmpDataQpScenario.rhsVec = (double*) malloc( sizeof(double) * rhsSize);
		//
		// set dimensions for symmetric matrix possibly adding diagonal elements 
		// 
		int totalNzQuadCost = tmpDataQpScenario.totalNzQuadCost;
		//
		// alloc symmetric data
		//
		int sizeSymmetricMatrix = totalNzQuadCost;
		sizeSymmetricMatrix += tmpDataQpScenario.totalDecisionVars;
		tmpDataQpScenario.totalNzSymQuadCost = 0;
		tmpDataQpScenario.idxRowNzSymQuadCost = (int*) malloc( sizeof(int) * sizeSymmetricMatrix);
		tmpDataQpScenario.idxColNzSymQuadCost = (int*) malloc( sizeof(int) * sizeSymmetricMatrix);
		tmpDataQpScenario.valNzSymQuadCost = (double*) malloc( sizeof(double) * sizeSymmetricMatrix);
		// 
		// init symmetric data
		// NOTE: tikhonov contribution is dealth with here
		// NOTE: it is added here
		// 
		double epsilon = tmpDataQpScenario.muTarget;
		double tikhonovMu = tmpDataQpScenario.tikhonovMu;
		int totalElementsDiagonal = 0;
		int i = 0;
		for (int line = 1; line <= tmpDataQpScenario.totalDecisionVars; ++line) {
			//
			// controls to identify dinagonal
			//
			int totalElementsThisLine = 0;
			bool foundDiagLine = false;
			//
			// iterate
			//
			while (i < totalNzQuadCost && nzRowsQuadPerturbationPerScenario_[globalIndexQuadPerturbation + i] == line) {
				// 
				// get data
				// 
				int row = nzRowsQuadPerturbationPerScenario_[globalIndexQuadPerturbation + i];
				int col = nzColsQuadPerturbationPerScenario_[globalIndexQuadPerturbation + i];
				double val = nzValsQuadPerturbationPerScenario_[globalIndexQuadPerturbation + i];
				//
				// test if should introduce diagonal element here
				//
				if (foundDiagLine == false && col > line) {
					tmpDataQpScenario.idxRowNzSymQuadCost[tmpDataQpScenario.totalNzSymQuadCost] = line;
					tmpDataQpScenario.idxColNzSymQuadCost[tmpDataQpScenario.totalNzSymQuadCost] = line;
					tmpDataQpScenario.valNzSymQuadCost[tmpDataQpScenario.totalNzSymQuadCost] = 0.0;
					tmpDataQpScenario.totalNzSymQuadCost += 1;
					totalElementsDiagonal += 1;
					foundDiagLine = true;
				}
				//
				// test if is upper triangle
				//
				if (col > row) {
					tmpDataQpScenario.idxRowNzSymQuadCost[tmpDataQpScenario.totalNzSymQuadCost] = row;
					tmpDataQpScenario.idxColNzSymQuadCost[tmpDataQpScenario.totalNzSymQuadCost] = col;
					tmpDataQpScenario.valNzSymQuadCost[tmpDataQpScenario.totalNzSymQuadCost] = val;
					tmpDataQpScenario.totalNzSymQuadCost += 1;
				} else if (col == row) {
					tmpDataQpScenario.idxRowNzSymQuadCost[tmpDataQpScenario.totalNzSymQuadCost] = row;
					tmpDataQpScenario.idxColNzSymQuadCost[tmpDataQpScenario.totalNzSymQuadCost] = col;
					tmpDataQpScenario.valNzSymQuadCost[tmpDataQpScenario.totalNzSymQuadCost] = val + (epsilon*tikhonovMu);
					tmpDataQpScenario.totalNzSymQuadCost += 1;
				}
				// 
				// count elements in diagonal
				// 
				if (col == row) {
					totalElementsDiagonal += 1;
					foundDiagLine = true;
				}
				//
				// counter
				//
				i += 1;
				totalElementsThisLine += 1;
			}
			//
			// in case of no elements this line
			//
			if (totalElementsThisLine == 0) {
				tmpDataQpScenario.idxRowNzSymQuadCost[tmpDataQpScenario.totalNzSymQuadCost] = line;
				tmpDataQpScenario.idxColNzSymQuadCost[tmpDataQpScenario.totalNzSymQuadCost] = line;
				tmpDataQpScenario.valNzSymQuadCost[tmpDataQpScenario.totalNzSymQuadCost] = epsilon*tikhonovMu;
				tmpDataQpScenario.totalNzSymQuadCost += 1;
				totalElementsDiagonal += 1;
			}
		}
		// 
		// check for error
		// 
		if (totalElementsDiagonal != tmpDataQpScenario.totalDecisionVars) {
			printf("totalElementsDiagonal: %d\n", totalElementsDiagonal);
			printf("totalDecisionVars: %d\n", tmpDataQpScenario.totalDecisionVars);
			printf("Do not has complete diagonal, scenario = %d.\n", s);
			exit(-1);
		}
		// 
		// set indexes
		// 
		globalIndexCostScenario += tmpDataQpScenario.totalDecisionVars;
		globalIndexQuadPerturbation += tmpDataQpScenario.totalNzQuadCost;
		globalIndexMatrixA += tmpDataQpScenario.totalNzMatrixA;
		//
		// append second stage data
		//
		data_->secondStageData.push_back(tmpDataQpScenario);
	}
	//
	// init buffer for solution and multipliers of second stage
	// 
	data_->optPrimalSolScen = (double*) malloc( sizeof(double) * totalDecisionVarsSecondStage);
	data_->optMultEqConstraintsScen = (double*) malloc( sizeof(double) * totalEqConstraintsSecondStage);
	//
	// init buffer for grad and val of second stage agregated
	//
	data_->lastFirstStageTrial = (double*) malloc( sizeof(double) * totalDecisionVarsFirstStage_);
	data_->lastValuesPerScenario = 10000000.0;
	data_->lastGradSecondStage = (double*) malloc( sizeof(double) * totalDecisionVarsFirstStage_);
	for (int i=0; i<totalDecisionVarsFirstStage_; ++i) {
		data_->lastFirstStageTrial[i] = -1.0;
		data_->lastGradSecondStage[i] = -1.0;
	}
	// 
	// init buffer for hessian values of second stage
	//
	if (data_->useSecondOrderInfoSecondStage == true) {
		int tmp = totalDecisionVarsFirstStage_;
		int totalPositions = tmp*(tmp + 1)/2;
		data_->totalPositionsHessSecondStage = totalPositions;
		//
		// alloc buffer
		//
		data_->idxRowHessSecondStage = (int*) malloc( sizeof(int) * totalPositions);
		data_->idxColHessSecondStage = (int*) malloc( sizeof(int) * totalPositions);
		data_->valsHessSecondStage = (double*) malloc( sizeof(double) * totalPositions);
		//
		// alloc buffer
		//
		int globalIndex = 0;
		for (int row=1; row<=totalDecisionVarsFirstStage_; ++row) {
			for (int col=row; col<=totalDecisionVarsFirstStage_; ++col) {
				data_->idxRowHessSecondStage[globalIndex] = row;
				data_->idxColHessSecondStage[globalIndex] = col;
				globalIndex += 1;
			}
		}
	}
	//
	// global counters to size buffers correctly
	//	
	int maxTotNewtonVars = -1;
	int maxTotalDecisionVars = -1;
	int maxTotalEqualityConstraints = -1;
	//
	// iterate scenarios to pre-compute newton matrices
	// 
	for (int s=0; s<totalScenarios_; ++s) {
		//
		// main data structure
		//
		struct QpProbData* dataScenIter = &(data_->secondStageData[s]);
		//
		// dimensions scenario
		//
		int totalDecisionVars = dataScenIter->totalDecisionVars;
		int totalEqualityConstraints = dataScenIter->totalEqualityConstraints;
		int totalNzQuadSymCost = dataScenIter->totalNzSymQuadCost;
		int totalNzMatrixA = dataScenIter->totalNzMatrixA;
		int totNewtonVars = totalDecisionVars + totalEqualityConstraints;
		//
		// max
		//
		if (totalDecisionVars > maxTotalDecisionVars) {
			maxTotalDecisionVars = totalDecisionVars;
		}
		if (totalEqualityConstraints > maxTotalEqualityConstraints) {
			maxTotalEqualityConstraints = totalEqualityConstraints;
		}
		//
		// alloc info to newton matrix iter
		//
		int totNzNewtonMatrix = totalNzQuadSymCost + totalNzMatrixA + totalEqualityConstraints;
		MKL_INT* rowIncrementNewtonMatrix = (MKL_INT*) malloc(sizeof(MKL_INT) * (totNewtonVars + 1));
		MKL_INT* colIdxNewtonMatrix = (MKL_INT*) malloc(sizeof(MKL_INT) * totNzNewtonMatrix);
		double* valsNewtonMatrix = (double*) malloc(sizeof(double) * totNzNewtonMatrix);
		int* positionsToChange = (int*) malloc(sizeof(int) * totalDecisionVars);
		//
		// fill data of lower triangular part of newton matrix to send to pardiso mkl
		//
		int globalIndexValNewtonMatrix = 0;
		int globalIndexSymQuadCost = 0;
		int globalIndexMatrixAt = 0;
		rowIncrementNewtonMatrix[0] = 1;
		//
		// iterate
		//
		for (int currentLineNewtonMatrix = 1; currentLineNewtonMatrix <= totalDecisionVars; ++currentLineNewtonMatrix) {
			//
			// init
			// 
			int totalNzElements = 0;
			//
			// quad cost
			//
			while (globalIndexSymQuadCost < totalNzQuadSymCost &&
					dataScenIter->idxRowNzSymQuadCost[globalIndexSymQuadCost] == currentLineNewtonMatrix) {
				int row = dataScenIter->idxRowNzSymQuadCost[globalIndexSymQuadCost];
				int col = dataScenIter->idxColNzSymQuadCost[globalIndexSymQuadCost];
				// test of diagonal element
				if (row == col) {
					positionsToChange[row-1] = globalIndexValNewtonMatrix;
				}
				// normal things
				colIdxNewtonMatrix[globalIndexValNewtonMatrix] = col;
				double tmp = dataScenIter->valNzSymQuadCost[globalIndexSymQuadCost];
				valsNewtonMatrix[globalIndexValNewtonMatrix] = 2*tmp;
				globalIndexSymQuadCost += 1;
				globalIndexValNewtonMatrix += 1;
				totalNzElements += 1;
			}
			//
			// matrix At
			//
			while (globalIndexMatrixAt < totalNzMatrixA &&
					dataScenIter->idxRowNzMatrixAt[globalIndexMatrixAt] == currentLineNewtonMatrix) {
				colIdxNewtonMatrix[globalIndexValNewtonMatrix] = 
					totalDecisionVars + dataScenIter->idxColNzMatrixAt[globalIndexMatrixAt];
				valsNewtonMatrix[globalIndexValNewtonMatrix] = dataScenIter->valNzMatrixAt[globalIndexMatrixAt];
				globalIndexMatrixAt += 1;
				globalIndexValNewtonMatrix += 1;
				totalNzElements += 1;
			}
			// 
			// row increment
			//
			int indexTmp = currentLineNewtonMatrix;
			rowIncrementNewtonMatrix[indexTmp] = totalNzElements + rowIncrementNewtonMatrix[indexTmp-1];
		}
		// 
		// max newton vars to alloc space for rhs and sol
		//
		if (maxTotNewtonVars < totNewtonVars) {
			maxTotNewtonVars = totNewtonVars;
		}
		// 
		// other lines are empty
		// 
		for (int i = totalDecisionVars+1; i <= totNewtonVars; ++i) {
			rowIncrementNewtonMatrix[i] = rowIncrementNewtonMatrix[i-1] + 1;
			colIdxNewtonMatrix[globalIndexValNewtonMatrix] = i;
			valsNewtonMatrix[globalIndexValNewtonMatrix] = 0.0;
			globalIndexValNewtonMatrix += 1;
		}
		//
		// append to pre-computed newton information
		//
		data_->sizesNewtonSystemsPerScenario.push_back(totNewtonVars);
		data_->rowIncrementsNewtonSystemsPerScenario.push_back(rowIncrementNewtonMatrix);
		data_->colsNewtonSystemsPerScenario.push_back(colIdxNewtonMatrix);
		data_->valsNewtonSystemsPerScenario.push_back(valsNewtonMatrix);
		data_->positionsOfDiagonalElementsToChangePerScenario.push_back(positionsToChange);
		
	}
	//
	// init space to store first order information
	//
	if (data_->useSecondOrderInfoSecondStage == true) {
		//
		// determine max size first stage vars
		//
		int maxSizePrimalVarsSecondStage = -1;
		for (int s=0; s<totalScenarios_; ++s) {
			int tmp = totalDecisionVarsPerScenario_[s];
			if (tmp > maxSizePrimalVarsSecondStage) {
				maxSizePrimalVarsSecondStage = tmp;
			}
		}
		//
		// alloc space for derivatives
		//
		int len = data_->totalDecisionVarsFirstStage;
		for (int i=0; i<len; ++i) {
			double* spaceToStoreDer = (double*) malloc( sizeof(double) * maxSizePrimalVarsSecondStage);
			data_->derivativesPrimalPerScenarioWrtPcoord.push_back(spaceToStoreDer);
		}
	}
	//
	// append references for reserved space
	// 
	data_->rhsForNewtonPerScenario = (double*) malloc( sizeof(double) * maxTotNewtonVars );
	data_->solForNewtonPerScenario = (double*) malloc( sizeof(double) * maxTotNewtonVars );
	//
	// alloc vectors for the refinement of the solution
	//
	data_->gradientOfLagrangianIter = (double*) malloc(sizeof(double) * maxTotNewtonVars);
	data_->gradLagrTrialPoint = (double*) malloc(sizeof(double) * maxTotNewtonVars);
	data_->newtonDirectionIter = (double*) malloc(sizeof(double) * maxTotNewtonVars);
	data_->primalSolIter = (double*) malloc(sizeof(double) * maxTotalDecisionVars);
	data_->primalSolIterDamp = (double*) malloc(sizeof(double) * maxTotalDecisionVars);
	data_->multiplierIter = (double*) malloc(sizeof(double) * maxTotalEqualityConstraints); 
	data_->multIterDamp = (double*) malloc(sizeof(double) * maxTotalEqualityConstraints);
	//
	// control vars
	// NOTE: this is used to start the iteration at the master problem
	//
	data_->levelIter = 1;
	data_->scenarioIter = -1;
}

void freeDataSolverUpperSmoothing(
		DataSolverUpperSmoothing* data_) {
	//
	// free data first stage data structure
	//
	free(data_->firstStageData.idxRowNzSymQuadCost);
	free(data_->firstStageData.idxColNzSymQuadCost);
	free(data_->firstStageData.valNzSymQuadCost);
	// 
	// iterate scenarios
	//
	int totalScenarios_ = data_->totalScenarios;
	for (int s=0; s<totalScenarios_; ++s) {
		//
		// free symmetric matrix data
		//
		free(data_->secondStageData[s].rhsVec);
		free(data_->secondStageData[s].idxRowNzSymQuadCost);
		free(data_->secondStageData[s].idxColNzSymQuadCost);
		free(data_->secondStageData[s].valNzSymQuadCost);
		//
		// free pre-computed newton matrices
		// 
		free(data_->rowIncrementsNewtonSystemsPerScenario[s]);
		free(data_->colsNewtonSystemsPerScenario[s]);
		free(data_->valsNewtonSystemsPerScenario[s]);
		free(data_->positionsOfDiagonalElementsToChangePerScenario[s]);
	}
	//
	// free rhs and sol
	//
	free(data_->rhsForNewtonPerScenario);
	free(data_->solForNewtonPerScenario);
	//
	// free buffer space
	//
	free(data_->lastFirstStageTrial);
	free(data_->lastGradSecondStage);
	free(data_->optPrimalSolScen);
	free(data_->optMultEqConstraintsScen);
	// 
	// free space used to store first order derivatives of solution mappings
	// 
	if (data_->useSecondOrderInfoSecondStage == true) {
		int len = data_->totalDecisionVarsFirstStage;
		for (int i=0; i<len; ++i) {
			free(data_->derivativesPrimalPerScenarioWrtPcoord[i]);
		}
	}
	//
	// free buffer for hessian of second stage
	//
	if (data_->useSecondOrderInfoSecondStage == true) {
		free(data_->idxRowHessSecondStage);
		free(data_->idxColHessSecondStage);
		free(data_->valsHessSecondStage);
	}
	//
	// free refinement vectors
	// 
	free(data_->gradientOfLagrangianIter);
	free(data_->gradLagrTrialPoint);
	free(data_->primalSolIter);
	free(data_->primalSolIterDamp);
	free(data_->multiplierIter);
	free(data_->newtonDirectionIter);
	free(data_->multIterDamp);
}

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
		double* optMulEqConsOut) {
	try {
		//
		// data
		//
		DataSolverUpperSmoothing data_;
		//	
		// build solver object
		//
		initDataSolverUpperSmoothing(
			&data_,
			useSecondOrderInfoMasterProblem,
			maxIterPerSubProblem,
			totalNewtonRefinementStepsAfterIpopt,
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
			nzValsQuadPerturbationPerScenario);
		//
		// solve
		//
		solveMuTargetSubproblemWithIpopt(
			-1,
			&data_,
			optSolOut,
			optValOut,
			optMulEqConsOut);
		//
		// free
		//
		freeDataSolverUpperSmoothing(&data_);
	} catch(...) {
		return -1;
	}
	//
	// return
	// 
	return 0;
}
