/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dualbound.c
 * @brief  unit tests for SCIPupdateLocalDualbound() in the presolving stages
 * @author Youssouf Emine
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

#define NVARS 30
#define EPS   1e-6

/** GLOBAL VARIABLES **/
static SCIP_Real presolbound;

/** execution method of presolver that passes a dual bound in original space to SCIP */
static
SCIP_DECL_PRESOLEXEC(presolExecDualbound)
{
   SCIP_CALL( SCIPupdateLocalDualbound(scip, presolbound) );

   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

static SCIP* scip = NULL;

/** creates a binary problem with two knapsack rows that presolving does not solve
 *
 *  The optimal value is 391 + offset for maximization and 323 + offset for minimization.
 */
static
void createProb(
   SCIP_OBJSENSE         objsense,           /**< objective sense */
   SCIP_Real             offset              /**< objective offset */
   )
{
   SCIP_VAR* vars[NVARS];
   SCIP_Real avals[NVARS];
   SCIP_Real bvals[NVARS];
   SCIP_CONS* cons;
   SCIP_Real lhs;
   SCIP_Real rhs;
   char name[SCIP_MAXSTRLEN];
   int i;

   SCIP_CALL( SCIPcreateProbBasic(scip, "dualbound") );
   SCIP_CALL( SCIPsetObjsense(scip, objsense) );
   SCIP_CALL( SCIPaddOrigObjoffset(scip, offset) );

   for( i = 0; i < NVARS; ++i )
   {
      avals[i] = 20 + (7 * i) % 31;
      bvals[i] = 18 + (5 * i) % 29;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &vars[i], name, 0.0, 1.0, avals[i] + (11 * i) % 13, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(scip, vars[i]) );
   }

   lhs = objsense == SCIP_OBJSENSE_MAXIMIZE ? -SCIPinfinity(scip) : 300.5;
   rhs = objsense == SCIP_OBJSENSE_MAXIMIZE ? 300.5 : SCIPinfinity(scip);
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "a", NVARS, vars, avals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   lhs = objsense == SCIP_OBJSENSE_MAXIMIZE ? -SCIPinfinity(scip) : 280.5;
   rhs = objsense == SCIP_OBJSENSE_MAXIMIZE ? 280.5 : SCIPinfinity(scip);
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "b", NVARS, vars, bvals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   for( i = 0; i < NVARS; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );
   }
}

/* TEST SUITE */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPincludePresolBasic(scip, NULL, "dualbound", "passes a dual bound to SCIP", 10000000, -1,
         SCIP_PRESOLTIMING_FAST, presolExecDualbound, NULL) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(dualbound, .init = setup, .fini = teardown);

/* TESTS */
Test(dualbound, presolving_max, .description = "tests dual bound passed in presolving for maximization")
{
   createProb(SCIP_OBJSENSE_MAXIMIZE, 0.0);
   presolbound = 1000.0;

   SCIP_CALL( SCIPpresolve(scip) );
   cr_expect_float_eq(SCIPgetDualbound(scip), 1000.0, EPS, "dual bound is %g", SCIPgetDualbound(scip));

   SCIP_CALL( SCIPsolve(scip) );
   cr_assert_eq(SCIPgetStatus(scip), SCIP_STATUS_OPTIMAL);
   cr_expect_float_eq(SCIPgetPrimalbound(scip), 391.0, EPS, "optimal value is %g", SCIPgetPrimalbound(scip));
}

Test(dualbound, presolving_min_offset, .description = "tests dual bound passed in presolving for minimization")
{
   createProb(SCIP_OBJSENSE_MINIMIZE, 500.0);
   presolbound = 400.0;

   SCIP_CALL( SCIPpresolve(scip) );
   cr_expect_float_eq(SCIPgetDualbound(scip), 400.0, EPS, "dual bound is %g", SCIPgetDualbound(scip));

   SCIP_CALL( SCIPsolve(scip) );
   cr_assert_eq(SCIPgetStatus(scip), SCIP_STATUS_OPTIMAL);
   cr_expect_float_eq(SCIPgetPrimalbound(scip), 823.0, EPS, "optimal value is %g", SCIPgetPrimalbound(scip));
}

Test(dualbound, presolved_max_offset, .description = "tests dual bound passed after presolving for maximization")
{
   createProb(SCIP_OBJSENSE_MAXIMIZE, 500.0);
   presolbound = 2000.0;

   SCIP_CALL( SCIPpresolve(scip) );
   cr_expect_float_eq(SCIPgetDualbound(scip), 2000.0, EPS, "dual bound is %g", SCIPgetDualbound(scip));

   SCIP_CALL( SCIPupdateLocalDualbound(scip, 1000.0) );
   cr_expect_float_eq(SCIPgetDualbound(scip), 1000.0, EPS, "dual bound is %g", SCIPgetDualbound(scip));

   SCIP_CALL( SCIPsolve(scip) );
   cr_assert_eq(SCIPgetStatus(scip), SCIP_STATUS_OPTIMAL);
   cr_expect_float_eq(SCIPgetPrimalbound(scip), 891.0, EPS, "optimal value is %g", SCIPgetPrimalbound(scip));
}
