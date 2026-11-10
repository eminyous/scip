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

/**@file   trysol.c
 * @brief  unit tests for the trysol primal heuristic
 * @author Youssouf Emine
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/heur_trysol.h"
#include "scip/nodesel_bfs.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_HEUR* heur;
static SCIP_VAR* x;

/* creates a problem with one variable 0 <= x <= 20 of objective coefficient 1 and only the trysol heuristic */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeHeurTrySol(scip) );
   SCIP_CALL( SCIPincludeNodeselBfs(scip) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );

   heur = SCIPfindHeur(scip, "trysol");
   cr_assert(heur != NULL);

   SCIP_CALL( SCIPcreateProbBasic(scip, "trysol") );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 20.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
}

/* releases variable, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

/** passes solutions x = 10 and x = 5 to the heuristic, runs it before presolving, and returns the committed value */
static
SCIP_Real passSols(
   SCIP_OBJSENSE         objsense,           /**< objective sense */
   SCIP_Bool             add                 /**< whether the solutions are passed to be added instead of tried */
   )
{
   SCIP_Real vals[2] = { 10.0, 5.0 };
   SCIP_SOL* sol;
   int i;

   SCIP_CALL( SCIPsetObjsense(scip, objsense) );
   SCIP_CALL( SCIPtransformProb(scip) );

   for( i = 0; i < 2; ++i )
   {
      SCIP_CALL( SCIPcreateOrigSol(scip, &sol, NULL) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, vals[i]) );

      if( add )
      {
         SCIP_CALL( SCIPheurPassSolAddSol(scip, heur, sol) );
      }
      else
      {
         SCIP_CALL( SCIPheurPassSolTrySol(scip, heur, sol) );
      }

      SCIP_CALL( SCIPfreeSol(scip, &sol) );
   }

   SCIP_CALL( SCIPpresolve(scip) );

   cr_assert_eq(SCIPgetNSols(scip), 1);
   sol = SCIPgetBestSol(scip);
   cr_assert_eq(SCIPsolGetHeur(sol), heur);

   return SCIPgetSolOrigObj(scip, sol);
}

TestSuite(trysol, .init = setup, .fini = teardown);

Test(trysol, trysolmax)
{
   cr_expect_eq(passSols(SCIP_OBJSENSE_MAXIMIZE, FALSE), 10.0);
}

Test(trysol, addsolmax)
{
   cr_expect_eq(passSols(SCIP_OBJSENSE_MAXIMIZE, TRUE), 10.0);
}

Test(trysol, trysolmin)
{
   cr_expect_eq(passSols(SCIP_OBJSENSE_MINIMIZE, FALSE), 5.0);
}

Test(trysol, addsolmin)
{
   cr_expect_eq(passSols(SCIP_OBJSENSE_MINIMIZE, TRUE), 5.0);
}
