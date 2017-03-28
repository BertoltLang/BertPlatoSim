#
# install.sh runs automatically this script. It can be done manually running (from the same folder than install.sh):
# $ python ./dependencies/installscripts/install_hdf5.py
#


import os,shutil,subprocess


# Specify the dependency package name

packageName = "hdf5-1.8.16"

# Specify build and install folders

currentWorkingDir = os.getcwd()
buildDir   = currentWorkingDir + "/dependencies/Downloads/"
installDir = currentWorkingDir + "/dependencies/Installs/" + packageName 

# Remove a possible older version, and create a fresh one

shutil.rmtree(installDir, ignore_errors=True)
os.mkdir(installDir)

# Print a banner

print("\n\n\n")
print("===============")
print("Installing HDF5")
print("===============")
print("\n")

# Build and install package

installProcedure = "cd {build};                                     \
					tar -xzvf {package}.tgz;                        \
					cd {package};                                   \
                    ./configure --prefix={install} --enable-cxx;    \
                    make;                                           \
                    make install".format(build=buildDir, package=packageName, install=installDir)

subprocess.call(installProcedure, shell=True)


# After installation in the install folder, remove the decompressed package folder in the build dir
# so that only the .tgz file remains in the Downloads folder.

shutil.rmtree(buildDir+packageName, ignore_errors=True)


