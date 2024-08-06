//SET BATCH MODE 
setBatchMode(true);   //true or false, true if you don't want to see the images, which is faster

//START MESSAGE
print("**** STARTING THE MACRO ****");

//INPUT/OUPUT folders
inDir=getDirectory("Choose the input folder"); 
outputDir=getDirectory("And the output folder");
myList=getFileList(inDir);  //an array

for (j = 0 ; j < myList.length ; j++ ){
	path=inDir+myList[j];   //path to each file
	open(path);
	FileName=File.nameWithoutExtension;
	ImageID=File.name;
	print("Processing "+ImageID);
	
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'"+ImageID+"', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.479071', 'nmsThresh':'0.3', 'outputType':'Label Image', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	selectWindow("Label Image");
	saveAs("Tiff", outputDir+ImageID);
	close("*");
}