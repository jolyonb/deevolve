import os

dirtoinis = "ginis"

DEevolexec = "./intDE "
for file in os.listdir(dirtoinis):
    if file.endswith(".gni"):
        print dirtoinis+"/"+file
	os.system(DEevolexec + dirtoinis+"/"+file)