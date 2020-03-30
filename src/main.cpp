/* Copyright 2020 Kristofer Björnson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/** @package TBTK2020_03_30
 *  @file main.cpp
 *  @brief New project
 *
 *  Empty template project.
 *
 *  @author Kristofer Björnson
 */

#include "TBTK/Model.h"
#include "TBTK/Range.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Streams.h"
#include "TBTK/TBTK.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Visualization/MatPlotLib/Plotter.h"

#include <complex>

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

int main(int argc, char **argv){
	//Initialize TBTK.
	Initialize();

	//Set the units.
	UnitHandler::setScales(
		{"1 rad", "1 C", "1 pcs", "1 J", "1 m", "1 K", "1 s"}
	);

	//Parameters.
	unsigned int N = 99;
	double a = 1e-9;
	Range x(-a/2, a/2, N+2);
	double hbar = UnitHandler::getConstantInNaturalUnits("hbar");
	double m_e = UnitHandler::getConstantInNaturalUnits("m_e");
	double t = hbar*hbar/(2*m_e*pow(x[1] - x[0], 2));

	//Set up the model.
	Model model;
	for(unsigned int n = 0; n < N; n++){
		model << HoppingAmplitude(2*t, {n}, {n});
		if(n + 1 < N)
			model << HoppingAmplitude(-t, {n}, {n+1});

		model << HoppingAmplitude(0.5*5000*x[n+1]*x[n+1], {n}, {n}) + HC;
	}
	model.construct();

	//Set up and run the solver.
	Solver::Diagonalizer solver;
	solver.setModel(model);
	solver.run();

	//Set up the property extractor.
	PropertyExtractor::Diagonalizer propertyExtractor(solver);

	//Extract eigenvalues and wave functions.
	Property::EigenValues eigenValues = propertyExtractor.getEigenValues();
	Property::WaveFunctions waveFunctions
		= propertyExtractor.calculateWaveFunctions(
			{{_a_}},
			{_a_}
		);

	//Plot the eigenvalues and wave functions.
	Plotter plotter;
	plotter.setTitle("Eigenvalues");
	plotter.setLabelX("Eigenvalue number");
	plotter.setLabelY("Energy");
	plotter.plot(eigenValues.getData());
	plotter.save("figures/EigenValues.png");
	for(unsigned int n = 0; n < N; n++){
		plotter.clear();
		plotter.plot({_a_}, n, waveFunctions);
		plotter.save("figures/WaveFuntion" + to_string(n) + ".png");
	}

	return 0;
}
