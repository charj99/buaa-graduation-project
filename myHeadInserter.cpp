#include "rose.h"
using namespace std;
using namespace SageInterface;
using namespace SageBuilder;

void myInsertHeader(SgProject* project) {
	SgGlobal* globalscope = getFirstGlobalScope(project);
	insertHeader("physis/physis.h", PreprocessingInfo::before, 
					false, globalscope);
}

int main (int argc, char** argv)
{
    // Build the AST used by ROSE
    SgProject* project = frontend(argc, argv);
    project->skipfinalCompileStep(true);
    // Insert your own manipulations of the AST here...		

	myInsertHeader(project);

	// Run internal consistency tests on AST
    AstTests::runAllTests(project);

    // Generate source code from AST and invoke your
    // desired backend compiler
    return backend(project);
}

