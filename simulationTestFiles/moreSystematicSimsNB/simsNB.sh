# I simulate 9 times using all combinations of {1,5,30} numMethods and numLocations
#!/bin/bash


# 1,5,30

for locs in 1 5 30;
	do for meths in 1 5 30;
		# simulate
		do outname=$(echo locs_"$locs"_meths_"$meths");
		/DATA/git/birp_cpp2/build/birp --task simulate --negativeBinomial --out /home/liam/ownCloud/institution/project_birp/simulationTestFiles/moreSystematicSimsNB/simulatedFiles/sims_"$outname" --timesOfChange 5 --gamma -0.5,0.5 --timepoints 1,2,3,4,5,6,7,8,9,10 --numLocations $locs --numMethods $meths;
		# infer
		/DATA/git/birp_cpp2/build/birp --task infer --negativeBinomial --counts /home/liam/ownCloud/institution/project_birp/simulationTestFiles/moreSystematicSimsNB/simulatedFiles/sims_"$outname"_simulatedCounts.txt --timesOfChange 5 --out /home/liam/ownCloud/institution/project_birp/simulationTestFiles/moreSystematicSimsNB/results/"$outname" --iterations 10000
	done;
done;

# /DATA/git/birp_cpp2/build/birp --task infer --negativeBinomial --counts /home/liam/ownCloud/institution/project_birp/simulationTestFiles/moreSystematicSimsNB/simulatedFiles/sims_locs_1_meths_1_simulatedCounts.txt --timesOfChange 5 --out /home/liam/ownCloud/institution/project_birp/simulationTestFiles/moreSystematicSimsNB/results/valgrind_sims_1_1_debug --iterations 10000


# /DATA/git/birp_cpp2/build/birp --task simulate --negativeBinomial --out /home/liam/ownCloud/institution/project_birp/simulationTestFiles/moreSystematicSimsNB/simulatedFiles/sims_locs_30_meths_30 --timesOfChange 5 --gamma -0.5,0.5 --timepoints 1,2,3,4,5,6,7,8,9,10 --numLocations 30 --numMethods 30 --sigmaMu 0.001

# /DATA/git/birp_cpp2/build/birp --task infer --negativeBinomial --counts /home/liam/ownCloud/institution/project_birp/simulationTestFiles/moreSystematicSimsNB/simulatedFiles/sims_locs_30_meths_30_simulatedCounts.txt --timesOfChange 5 --out /home/liam/ownCloud/institution/project_birp/simulationTestFiles/moreSystematicSimsNB/results/locs_30_meths_30 --iterations 10000
