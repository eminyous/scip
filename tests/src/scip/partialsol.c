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

/**@file   partialsol.c
 * @brief  unit test for the partial solution getters in scip_sol.c
 * @author Youssouf Emine
 */

#include "scip/scip.h"
#include "include/scip_test.h"

/** GLOBAL VARIABLES **/
static SCIP* scip = NULL;
static SCIP_VAR* var = NULL;

/* TEST SUITE */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "partialsol") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &var, "x", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, var) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &var) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(partialsol, .init = setup, .fini = teardown);

Test(partialsol, count, .description = "tests that added partial solutions are counted and returned")
{
   SCIP_SOL* sol;
   SCIP_Bool stored;

   cr_expect_eq(SCIPgetNPartialSols(scip), 0);

   SCIP_CALL( SCIPcreatePartialSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, var, 1.0) );
   SCIP_CALL( SCIPaddSolFree(scip, &sol, &stored) );

   cr_assert(stored);
   cr_assert_eq(SCIPgetNPartialSols(scip), 1);
   cr_assert_not_null(SCIPgetPartialSols(scip));
   cr_expect(SCIPsolIsPartial(SCIPgetPartialSols(scip)[0]));
   cr_expect_eq(SCIPgetSolVal(scip, SCIPgetPartialSols(scip)[0], var), 1.0);
}
