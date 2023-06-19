//option 1
setBatchMode("hide");
setThreshold(100,100000000000000); //put your thresholds here
for (i = 1; i <= nSlices; i++) {
    setSlice(i);   
	run("Create Selection");
	run("Set...", "value=-1000");
	run("Select None");
}
resetThreshold;
setBatchMode("exit and display");