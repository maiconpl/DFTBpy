import os
import codecs
import re
from distutils.core import setup
# https://docs.python.org/2/distutils/setupscript.html#installing-package-data

# --BEGIN -- to be edited by the user (installer):
compiler = "gcc"
# --END

base_dir = os.path.abspath(os.path.dirname(__file__))

def compile_library():
    path = base_dir + "/src/lib"
    list_of_code = ["/build_dgamma_matrix.c", "/build_gamma_matrix.c", "/build_P_matrix.c", "/build_scc_matrix.c", "/polint.c"]
    list_of_objects = ["/build_dgamma_matrix.so", "/build_gamma_matrix.so", "/build_P_matrix.so", "/build_scc_matrix.so", "/polint.so"]
    # The user has to choose one optimization flag for the gcc compiler:
    compiler_flags = " -shared -fPIC -g -o "
#    compiler_flags = " -O3 -shared -fPIC -g -o "
#    compiler_flags = "-O2 -ffast-math -fexpensive-optimizations -shared -fPIC -g -o "
    for i in xrange(len(list_of_code)):
         all_compiler_flags = compiler_flags + path + list_of_objects[i] + " " + path + list_of_code[i]
         os.system("rm " + path + list_of_objects[i])
         os.system(compiler + all_compiler_flags)

compile_library()

def find_version(*file_paths):
    # Open in Latin-1 so that we avoid encoding errors.
    # Use codecs.open for Python 2 compatibility
    with codecs.open(os.path.join(base_dir, *file_paths), 'r', 'latin1') as f:
        version_file = f.read()

    # The version line must have the form
    # __version__ = 'ver'
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

compile_library()

setup(
    name = "DFTBpy",
    version = find_version("src", "__init__.py"),
    author = "Maicon Pierre Lourenco (MPL)",
    author_email = "maiconpl01@gmail.com",
    description = "It performs DFTB/SCC-DFTB calculation for molecules.",
    license = "GPL",
    keywords = "DFTB, SCC-DFTB, geometry optimization",
    packages = [
       "src",
       "src.model",
       "src.view",
       "src.util",
       "src.lib"
    ],
    package_data = {"src": ["lib/*c", "lib/*so"]}
)

def create_executable():
    dirname = os.listdir(base_dir + "/build")
    os.system("cp " + base_dir + "/build/" + dirname[0] + "/src/__init__.py " + base_dir + "/build/" + dirname[0] + "/src/dftbpy")
    os.system("chmod +x " + base_dir + "/build/" + dirname[0] + "/src/dftbpy")

create_executable()
