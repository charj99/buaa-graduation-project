#include "rose.h"
#include <assert.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <boost/lexical_cast.hpp>
#include <getopt.h>
#include <cstdio>
using namespace std;
using namespace SageInterface;
using namespace SageBuilder;
#define DEBUG 1
#define dbg if (DEBUG) cerr
#define GAP "--------------------\n"

#define SgForStatement_BODY_INDEX 3
#define FOREACH(type, list, var) \
    for (type::iterator var = list.begin(); var != list.end(); var++) 

typedef Rose_STL_Container<SgNode*> SgNodeList;
typedef set<SgNode*> SgNodeSet;
typedef pair<SgNode*, SgNodeList> ArrayInfo;
typedef vector<ArrayInfo> ArrayInfoList;
typedef map<SgName, SgType*> VarInfoMap;
typedef map<SgName, SgNodeList> ArrayRangeMap;
typedef pair<SgNode*, SgNode*> myRange;
typedef vector<pair<SgName, myRange> > RangeList;
typedef vector<SgExpression*> SgExpressionList;
typedef map<SgName, SgTypedefType*> SgTypedefTypeMap;
typedef map<SgName, SgName> CircularArrMap; 
typedef vector<string> StrList;
typedef vector<SgName> FieldList;
typedef map<SgName, FieldList> FieldRangeMap;
typedef map<SgType*, FieldRangeMap> StructFieldMap; 
// typedef vector<SgVarRefExp*> SgVarRefList;

/* structures to record stencil information */
struct DomainInfo {
    RangeList d_range;          /* range of the domain */
    SgNodeList d_exprs;         /* expressions in a domain */
    VarInfoMap d_vars;          /* used variables in a domain */
    ArrayRangeMap d_arrays;     /* used arrays in a domain */
};

typedef vector<DomainInfo> DomainInfoList;

struct IterateInfo {
    myRange it_times;           /* iterator times of the stencil */
    bool it_circular;           /* wether using a cirular array in a stencil */
    CircularArrMap it_cirarrs;  /* circular array pairs in a stencil */
};

struct StencilInfo {
    SgNodeList st_nodes;        /* for statements which stand for domains in a stencil */
    int st_dim;                 /* dimension of the stencil */
    SgNodeList st_dimMaxSize;   /* max range of each dimension */
    ArrayRangeMap st_arrays;    /* used arrays for a stencil */
    enum ST_TYPE {
        ST_TYPE_UNKNOWN = -1,
        ST_TYPE_UNSET,
        ST_TYPE_DOUBLE,
        ST_TYPE_FLOAT,
    } st_type;                  /* stencil type */
    DomainInfoList st_domains;  /* different kernels in a stencil */
    IterateInfo st_iter;        /* infomation about iterate times */
};

typedef vector<StencilInfo> StencilInfoList;
typedef map<SgNode*, StencilInfo> ForStencilInfoMap;

namespace physis {
// ##############################
//          find stencils
// ##############################
/* query ALL for statement */
void getForStmtList(SgNodeList &forStmtList, 
	SgProject* project) {
	
	dbg << "\n" << GAP << "ALL for statements\n" << GAP;	
	
	forStmtList = NodeQuery::querySubTree(project, V_SgForStatement);
    FOREACH(SgNodeList, forStmtList, iter) {
        SgForStatement* forStmt = isSgForStatement(*iter);
		
		dbg << forStmt->unparseToString() + "\n\n";
	}
    dbg << "# TOTALLY " << forStmtList.size() << " FOR NESTS\n\n";
}

bool isOutestFor(const SgNode* node) {
	SgNode* pNode = node->get_parent();	
	if (isSgForStatement(pNode)) return false;
	if (isSgBasicBlock(pNode) 
		&& isSgForStatement(pNode->get_parent())) return false;
	return true;
}

/* keep OUTEST for statements only */
void getOutestForList(SgNodeList &outestForList, 
	SgNodeList &forStmtList) {
	
	dbg << GAP << "OUTEST for statements\n" << GAP;	
    
    FOREACH(SgNodeList, forStmtList, iter) {
		SgNode* node = *iter;
		if (!isOutestFor(node)) continue;
		outestForList.push_back(node);	
		
		dbg << node->unparseToString() + "\n\n";
	}
    dbg << "# TOTALLY " << outestForList.size() << " FOR NESTS\n\n";
}

/* get loop body for for statement */
void getLoopBody(SgNodeList &exprList, SgNode* node) {
    SgForStatement* forNode = isSgForStatement(node);
    if (!forNode) return;
    SgNode* bodyNode = node->get_traversalSuccessorByIndex(
            SgForStatement_BODY_INDEX);
    SgBasicBlock* bbNode;
    if ((bbNode = isSgBasicBlock(bodyNode))) 
        exprList = bbNode->get_traversalSuccessorContainer();
    else exprList.push_back(bodyNode);
}

/*
 * check wether the for statement have a inner for statement
 * if so, return the for statement
 * else return NULL
 */
SgNode* haveInnerFor(SgNode* node) {
    SgNodeList exprList;
    getLoopBody(exprList, node);
    SgNode* forNode;
    FOREACH(SgNodeList, exprList, iter) {
        if ((forNode = isSgForStatement(*iter))) return forNode;
    }
    return NULL;
}

/* get the innest for loop */
SgNode* getInnestFor(SgNode* node) {
    SgNode* outerFor = node;
    SgNode* innerFor;
    while ((innerFor = haveInnerFor(outerFor))) outerFor = innerFor;
    return outerFor;
}

/* 
 * get information for the array visit, 
 * eg. a[i][j][k]
 */
void getArrayInfo4Arr(ArrayInfo &info, SgNode* node) {
    SgPntrArrRefExp* arrRef = isSgPntrArrRefExp(node);
    if (!arrRef) {
        info.first = node;
        return;
    }
    
    // dbg << "[array reference] " << arrRef->unparseToString() << "\n";
    // dbg << "[lhs of arrRef] " << arrRef->get_lhs_operand()->unparseToString() << "\n";
    // dbg << "[rhs of arrRef] " << arrRef->get_rhs_operand()->unparseToString() << "\n"; 
    
    info.second.push_back(arrRef->get_rhs_operand());
    getArrayInfo4Arr(info, arrRef->get_lhs_operand());
}

/* 
 * get information for the visited array in expression, 
 * eg. stencil: a[i][j][k] = c * a[i + 1][j][k],
 * get a[i][j][k] and a[i + 1][j][k]
 */
void getArrayInfo4Expr(ArrayInfoList &infoList, SgNode* node) {
    SgPntrArrRefExp* arrRef = isSgPntrArrRefExp(node);
    if (arrRef) {
        ArrayInfo info;
        getArrayInfo4Arr(info, node);
        infoList.push_back(info);
        return;
    }
    
    /* a[i][j][k] = b[i][j][k] = ... */
    if (isSgAssignOp(node)) return;
    
    SgNodeList children = node->get_traversalSuccessorContainer();
    FOREACH(SgNodeList, children, iter) {
        getArrayInfo4Expr(infoList, *iter);
    }
}

void printArrayInfoList(ArrayInfoList &infoList, bool isLhs = true) {
    if (infoList.empty()) return;

    if (isLhs) {
        dbg << "[lhs array info]\n";
    } else dbg << "[rhs array info]\n";
    

    FOREACH(ArrayInfoList, infoList, iter) {
        dbg << "\t[array identifier] " << (*iter).first->unparseToString() << "\n";
        SgNodeList indexs = (*iter).second;
        for (int i = (int)indexs.size() - 1; i >=0; i--) {
            dbg << "\t[array index] " << indexs[i]->unparseToString() << "\n";
        }
    }
}

/* judge wether a and b refer to the same variable */
bool sameValRef(SgNode* aNode, SgNode* bNode) {
    SgVarRefExp *a = isSgVarRefExp(aNode);
    SgVarRefExp *b = isSgVarRefExp(bNode);
    return a && b && a->get_symbol()->get_name() == b->get_symbol()->get_name(); 
}

/* judge wether b[]...[] is a neighbor of a[]...[] */
bool isNeighbor(ArrayInfo &a, ArrayInfo &b) {
    SgNodeList aIndexList = a.second;
    SgNodeList bIndexList = b.second;

    int dim = aIndexList.size();
    if ((int)bIndexList.size() != dim) {
        dbg << "\t[some rhs does not have same dimension with lhs]\n";
        return false;
    }

    SgBinaryOp* opNode;
    SgExpression *lhs, *rhs;
    bool result = true;
    for (int i = dim - 1; i >= 0; i--) {
        /* a[i] = ... b[i] ... */
        if (isSgVarRefExp(bIndexList[i])) { 
            result = sameValRef(aIndexList[i], bIndexList[i]);
            if (!result)
                dbg << "\t[rhs does not refer to the same variable with lhs]\n";
            return result;
        }
        else if ((opNode = isSgAddOp(bIndexList[i])) 
                || (opNode = isSgSubtractOp(bIndexList[i]))) {
            lhs = opNode->get_lhs_operand();
            rhs = opNode->get_rhs_operand();
            if (isSgCastExp(rhs))
                rhs = isSgCastExp(rhs)->get_operand();
            result = sameValRef(aIndexList[i], (SgNode* )lhs) && isSgIntVal((SgNode *)rhs);
                /* || sameValRef(bIndexList[i], (SgNode* )rhs) && isSgIntVal((SgNode *)lhs) */
            if (!result) 
                dbg << "\t[the index of rhs is not a neighbor of lhs]\n";
            return result;
        } else {
            dbg << "\t[unknown fault in isNeighbor()]\n";
            return false;
        }
    }
    return true; /* useless */
}

/* judge wether the expression is a stenicl assign operation */
bool isStencilAssign(SgNode *node) {
    SgAssignOp* assignNode = isSgAssignOp(node);
    if (!assignNode) {
        dbg << "\t[not an assign expression]\n";
        return false;
    }
    
    dbg << "[assign expression]\n\t" << assignNode->unparseToString() << "\n";

    /* get array reference infomation for both operands */
    SgExpression *lhs = assignNode->get_lhs_operand();
    SgExpression *rhs = assignNode->get_rhs_operand();
    ArrayInfoList lhInfoList, rhInfoList;
    getArrayInfo4Expr(lhInfoList, (SgNode* )lhs);
    getArrayInfo4Expr(rhInfoList, (SgNode* )rhs);
    
    /* for debug */
    printArrayInfoList(lhInfoList);
    printArrayInfoList(rhInfoList, false);

    /* 
     * judge wehter all arrRefs in rhs are neighbors of arrRef in lhs
     * the expression must have format: 
     * a[i][j][k] = f(x[i + a][j + b][k + c]) 
     */
    if (lhInfoList.size() != 1 || rhInfoList.empty()) {
        dbg << "\t[something wrong with lhs or rhs of the assign expression]\n";
        return false;
    }
    ArrayInfo center = lhInfoList[0];
    /*
    for (int i = center.second.size() - 1; i >= 0; i--) {
        SgVarRefExp* var = isSgVarRefExp(center.second[i]);
        if (!var) return false;
    }
    */
    for (int i = 0; i < (int)rhInfoList.size(); i++) {
        if (!isNeighbor(center, rhInfoList[i])) {
            return false;
        }
    }

    return true;
}

/* judge wether the expression is a stencil expression */
bool isStencilExpr(SgNode* node) {
    /* assign expression */
    SgExprStatement* exprStmt = isSgExprStatement(node);
    if (!exprStmt) {
        dbg << "\t[not a expression statement]\n";
        return false;
    }
    SgExpression* expr = exprStmt->get_expression(); 
    return isStencilAssign((SgNode *)expr) /* || TODO: other type of expressions */;
}

/* judge wether the for statement is a stencil */
bool isStencilForNest(SgNode* node) {
    /* 1: inner for statement must access neighbor gird */
    
    /* find innest for statement */
    SgNode* innestFor = getInnestFor(node);
    
    dbg << "[innest for]\n\t" << innestFor->unparseToString() << "\n";
    
    /* get expressions of the innest for statement */
    SgNodeList exprList;
    getLoopBody(exprList, innestFor);
    
    /* all expressions must(TODO: or may?) be stencil expression */
    FOREACH(SgNodeList, exprList, iter) {
        dbg << "[innest for body]\n\t" << (*iter)->unparseToString() << "\n";
        if (!isStencilExpr(*iter)) {
            return false;
        }
    }   

    /* 2: TODO */
    return true;
}

/* keep STENCILS only */
void getStencilNodeList(SgNodeList &stencilNodeList, 
        SgNodeList &outestForList) {
	
	dbg << GAP << "finding STNECILS\n" << GAP;	
    
    int count = 0;
    FOREACH(SgNodeList, outestForList, iter) {
        SgNode* node = *iter;
        dbg << "[for nest " << ++count << "]\n"; 
        if (!isStencilForNest(node)) { 
            dbg << "\n";
            continue;
        }
        dbg << "[STENCIL]\n";
        stencilNodeList.push_back(node);
        
        dbg << "\n";
    }
    dbg << "# TOTALLY " << outestForList.size() << " FOR NESTS\n";
    dbg << "# INCLUDING " << stencilNodeList.size() << " STENCILS\n\n"; 
}
// ##############################
//      get stencil info
// ##############################
/* print stencil infomation */
void printStencilInfo(StencilInfo* info, int stencilId) {  
    dbg << "[stencil " << stencilId << "]\n";
    
    dbg << "[location]\n";
    Sg_File_Info* fileInfo = info->st_nodes[0]->get_file_info();
    dbg << "\t" << fileInfo->get_filenameString();
    dbg << ": line " << fileInfo->get_line() 
        << ", col " << fileInfo->get_col() << "\n";

    dbg << "[st_nodes]\n";
    FOREACH(SgNodeList, info->st_nodes, iter) 
        dbg << "\t" << (*iter)->unparseToString() << "\n";
    dbg << "\n";

    dbg << "[st_dim]\n\t";
    dbg << info->st_dim << "D\n";

    dbg << "[st_dimMaxSize]\n";
    dbg << "\t(";
    FOREACH(SgNodeList, info->st_dimMaxSize, iter) {
        dbg << (*iter)->unparseToString() << ", ";
    }
    dbg << ")\n";

    dbg << "[st_type]\n\t";
    switch (info->st_type) {
        case StencilInfo::ST_TYPE_UNSET:
            dbg << "unset\n";
            break;
        case StencilInfo::ST_TYPE_DOUBLE:
            dbg << "double\n";
            break;
        case StencilInfo::ST_TYPE_FLOAT:
            dbg << "float\n";
            break;
        case StencilInfo::ST_TYPE_UNKNOWN:
        default:
            dbg << "unknown\n";
            break;
    }

    /* this is for array copy */
    dbg << "[st_arrays]\n";
    FOREACH(ArrayRangeMap, info->st_arrays, iter) {
        dbg << "\t" << iter->first.getString() << ": (";
        SgNodeList range = iter->second;
        for (int i = 0; i < (int) range.size(); i++) {
            dbg << range[i]->unparseToString() << ", "; 
        }
        dbg << ")\n";
    }
    
    int domainId = 0;
    dbg << "[st_domains]\n";
    FOREACH(DomainInfoList, info->st_domains, it) {
        dbg << "\t[domain " << ++domainId << "]\n";
	    dbg << "\t[d_range]\n";
	    FOREACH(RangeList, it->d_range, iter) {
	        dbg << "\t\t" << (*iter).first.getString() << ": ";
	        myRange range = (*iter).second;
	        dbg << "[" << (range.first)->unparseToString() 
	            << ", " << (range.second)->unparseToString() << ")\n";
	    }
	    dbg << "\t[d_exprs]\n";
	    FOREACH(SgNodeList, it->d_exprs, iter) {
	        dbg << "\t\t" << (*iter)->unparseToString() << "\n";
	    }
        dbg << "\t[d_vars]\n\t\t";
        FOREACH(VarInfoMap, it->d_vars, iter) {
            dbg << iter->first.getString() << ", ";
        }
        dbg << "\n";

	    dbg << "\t[d_arrays]\n";
	    FOREACH(ArrayRangeMap, it->d_arrays, iter) {
	        dbg << "\t\t" << iter->first.getString() << ": (";
	        SgNodeList range = iter->second;
	        for (int i = 0; i < (int) range.size(); i++) {
	            dbg << range[i]->unparseToString() << ", "; 
	        }
	        dbg << ")\n";
	    }
    }
    
    dbg << "[st_iter]\n";
    dbg << "\t[it_times]\n";
    myRange it_times = info->st_iter.it_times;
    dbg << "\t\t[" << (it_times.first)->unparseToString() 
            << ", " << (it_times.second)->unparseToString() << ")\n";
    dbg << "\t[it_circular]\n";
    if (info->st_iter.it_circular) { 
        dbg << "\t\ttrue\n";
    } else     
        dbg << "\t\tfalse\n";
    dbg << "\t[it_cirarrs]\n";
    CircularArrMap tmpMap = info->st_iter.it_cirarrs;
    FOREACH(CircularArrMap, tmpMap, iter) {
        dbg << "\t\t(" << iter->first.getString() 
            << ", " << iter->second.getString() << ")\n";
        tmpMap.erase(iter->second);
    }
}

SgNode* filterOutCast4Int(SgNode* node) {
    SgNode* result = deepCopyNode(node);    
    SgNodeList intList = NodeQuery::querySubTree(result, V_SgIntVal);
    FOREACH(SgNodeList, intList, iter) {
        SgExpression* curExp = isSgExpression(*iter);
        while (curExp && isSgCastExp(curExp->get_parent())) {
            curExp = isSgExpression(curExp->get_parent()); 
        }
        replaceExpression(curExp, isSgExpression(*iter));
    }
    return result;
}

/* 
 * get loopRange [a, b) 
 * in for statement for (int i = 1; i < MAXN; i++)
 * loopRangeList.push_back(make_pair(i, MAXN))
 */
void getLoopRange(RangeList &loopRangeList, SgNode *node) {
    SgForStatement* forNode = isSgForStatement(node);
    if (!forNode) return; /* something wrong */
    SgStatement* forInit = forNode->get_for_init_stmt();
    SgStatement* forTest = forNode->get_test();
    
    /* get identifier */
    SgNodeList idList = NodeQuery::querySubTree(forInit, V_SgInitializedName);
    SgInitializedName* initName = isSgInitializedName(idList[0]);
    SgName id = initName->get_name();
    
    /* get range start */
    SgNodeList stList = NodeQuery::querySubTree(initName, V_SgAssignInitializer);
    SgAssignInitializer* assignInit = isSgAssignInitializer(stList[0]);
    SgNode* st = assignInit->get_operand_i();
    if (isSgCastExp(st)) {
        st = isSgCastExp(st)->get_operand();
    }

    /* get range end */
    SgNodeList endList = NodeQuery::querySubTree(forTest, V_SgLessThanOp);
    SgLessThanOp* lessThan = isSgLessThanOp(endList[0]);
    if (!lessThan) { 
        /* TODO: less than or equal op */
        return;
    }
    SgNode* end = lessThan->get_rhs_operand();
    myRange range = make_pair(st, filterOutCast4Int(end));
    loopRangeList.push_back(make_pair(id, range));
}

/* get name in string for arrow exp, e.g. for fdtd->ez, return "fdtd->ez" */
string getName4ArrowExp(SgArrowExp *node) {
    if (!node) return "";
    SgVarRefExp* lhs = isSgVarRefExp(node->get_lhs_operand());
    SgVarRefExp* rhs = isSgVarRefExp(node->get_rhs_operand());
    string lName = lhs->get_symbol()->get_name();
    string rName = rhs->get_symbol()->get_name();
    return lName + "->" + rName;
}

/* 
 * get information for the visited variable in expression, 
 * eg. stenci: a[i][j][k] = c * a[i + 1][j][k],
 * get c
 */
void getVarInfo(VarInfoMap &infoMap, SgNode* node) {
    if (isSgPntrArrRefExp(node)) return;
    
    SgArrowExp* arrowExp = isSgArrowExp(node);
    if (arrowExp) {
        string varName = getName4ArrowExp(arrowExp);
        infoMap[varName] = arrowExp->get_type();
        return;
    }

    SgVarRefExp* varRef = isSgVarRefExp(node);
    if (varRef) {
        SgSymbol* symbol = varRef->get_symbol();
        infoMap[symbol->get_name()] = symbol->get_type();
        return;
    }
    
    SgNodeList children = node->get_traversalSuccessorContainer();
    FOREACH(SgNodeList, children, iter) {
        getVarInfo(infoMap, *iter);
    }
}

/* 
 * judge wether the stencil uses circular array
 * TODO: now only imply a naive judgement, need to be modified 
 */
/*
bool useCircularArray(SgNode* innestFor, SgNode* node) {
    SgNode* curNode = innestFor;
    while (curNode != node) {
        curNode = curNode->get_parent();
        if (isSgBasicBlock(curNode) && curNode->get_numberOfTraversalSuccessors() != 1) {
            return true;
        }
    }
    return false;
}*/

/* get alloc-site */
void getAllocSite(SgNodeSet &allocSites, SgProject* project) {
    SgNodeList funcCalls = NodeQuery::querySubTree(project, V_SgFunctionCallExp);
    FOREACH(SgNodeList, funcCalls, iter) {
        SgFunctionCallExp* funcCallExp = isSgFunctionCallExp(*iter);
        SgFunctionRefExp* funcRefExp 
            = isSgFunctionRefExp(funcCallExp->get_function());
        if (!funcRefExp) {
            /* something strange happend */
            continue;
        }
        SgFunctionSymbol* funcSymbol = funcRefExp->get_symbol();
        SgName funcName = funcSymbol->get_name();
        if (funcName == SgName("calloc")) 
            allocSites.insert(funcCallExp);
    }

    /* DEBUG: print alloc-sites */
    /*
    cout << GAP << "ALLOC-sites\n" << GAP;
    FOREACH(SgNodeSet, allocSites, iter) {
        cout << (*iter)->unparseToString() << "\n";
    }
    cout << "\n";
    */
}

/* judge wether the type of the initialized name is struct */
bool isStructType(SgType* type) {

    if (!isSgClassType(type))
        return false;

    SgClassType* classType = isSgClassType(type);
    if (!isSgClassDeclaration(classType->get_declaration())) 
        return false;

    SgClassDeclaration* classDecl 
        = isSgClassDeclaration(classType->get_declaration());
    if (classDecl->get_class_type() != SgClassDeclaration::e_struct) 
        return false;

    return true;
}

/* get the shape of dynamically allocated arrays in struct */
void getShape4AllocArrInStruct(SgExpressionList &expList,
        SgArrowExp* arrowExp, StructFieldMap &structFieldMap) {
    SgExpression* lhs = arrowExp->get_lhs_operand();
    SgExpression* rhs = arrowExp->get_rhs_operand();
    
    /* deal with lhs */
    if (!isSgVarRefExp(lhs)) return;
    SgInitializedName* initName 
        = isSgVarRefExp(lhs)->get_symbol()->get_declaration(); 
    SgType* type = initName->get_type();
    if (isSgPointerType(type)) 
        type = isSgPointerType(type)->get_base_type();
    if (!isStructType(type)) return;

    /* deal with rhs */
    if (!isSgVarRefExp(rhs)) return;
    SgSymbol* symbol = isSgVarRefExp(rhs)->get_symbol();
    SgName field = symbol->get_name();

    StructFieldMap::iterator it1 = structFieldMap.find(type);
    if (it1 == structFieldMap.end()) return;
    FieldRangeMap::iterator it2 = it1->second.find(field);
    if (it2 == it1->second.end()) return;

    FOREACH(FieldList, it2->second, iter) {
        expList.push_back(buildArrowExp(
                    copyExpression(lhs), 
                    buildVarRefExp(*iter, symbol->get_scope()))); 
    }
}

/* get the shape of the array */
SgNodeList getShape(SgNode* node, StructFieldMap &structFieldMap) {
    SgNodeList result;
	SgExpressionList expList;
    
    /* get shape of an array, e.g. a[i] */
    SgVarRefExp* varRef = isSgVarRefExp(node);
    if (varRef) {
	    SgInitializedName* initName = varRef->get_symbol()->get_declaration();
	    SgAssignInitializer* initer = isSgAssignInitializer(initName->get_initializer());
	    if (initer) {
	        /* declared and defined */
	        SgType* rhsType = initer->get_operand_i()->get_type();
            if (isSgArrayType(rhsType)) {
                expList = get_C_array_dimensions(*isSgArrayType(rhsType));
            } else { /* DOING */
                SgExpression* exp = initer->get_operand();
                while (isSgCastExp(exp))
                    exp = isSgCastExp(exp)->get_operand();
                getShape4AllocArrInStruct(expList, 
                        isSgArrowExp(exp), structFieldMap);
            }
	    } else {
	        /* declared only */
            if (isSgArrayType(initName->get_type())) {
	            expList = get_C_array_dimensions(*isSgArrayType(initName->get_type()));
            } /* else TODO: */
	    }
	    result.insert(result.end(), expList.begin(), expList.end());
	    return result;
    }
    
    /* get shape of an array of a struct, e.g. struct->a[i] */
    SgArrowExp* arrowExp = isSgArrowExp(node);
    if (arrowExp) {
        getShape4AllocArrInStruct(expList, 
                arrowExp, structFieldMap);
	    result.insert(result.end(), expList.begin(), expList.end());
        return result;
    }

    return result; // useless
}

/* get circular arrays */
void getCircularArrays(CircularArrMap& cirArrMap, SgForStatement* forStmt) {
    SgNodeList body;
    getLoopBody(body, forStmt);
    if (body.empty()) return;
   
    /* build transMap */
    CircularArrMap transMap;
    FOREACH(SgNodeList, body, iter) {
        if (isSgExprStatement(*iter)) {
            SgExpression* expr = isSgExprStatement(*iter)->get_expression();
            if (isSgAssignOp(expr)) {
                /* assign op, e.g. a = b; */
                SgAssignOp* assignOp = isSgAssignOp(expr);
                
                if (!isSgVarRefExp(assignOp->get_lhs_operand()))
                    continue;

                SgName lName = isSgVarRefExp(assignOp->get_lhs_operand())
                    ->get_symbol()->get_name();
                SgExpression* rExpr = assignOp->get_rhs_operand();
                if (isSgCastExp(rExpr))
                    rExpr = isSgCastExp(rExpr)->get_operand();
               
                if (!isSgVarRefExp(rExpr))
                    continue;

                SgName rName = isSgVarRefExp(rExpr)
                    ->get_symbol()->get_name();
                CircularArrMap::iterator it = transMap.find(rName);
                if (it != transMap.end()) 
                    transMap[lName] = it->second;
                else transMap[lName] = rName;
            }
        }
        if (isSgVariableDeclaration(*iter)) {
            /* variable declaration, e.g. void* a = b; */
            SgInitializedNamePtrList initNameList 
                = isSgVariableDeclaration(*iter)->get_variables();
            FOREACH(SgInitializedNamePtrList, initNameList, it) {
                SgName lName = (*it)->get_name();
                if (isSgAssignInitializer((*it)->get_initializer())) {
                    SgExpression* expr 
                        = isSgAssignInitializer((*it)->get_initializer())
                        ->get_operand();
                    if (isSgCastExp(expr)) 
                        expr = isSgCastExp(expr)->get_operand();
                
                    if (!isSgVarRefExp(expr))
                        continue;

                    SgName rName = isSgVarRefExp(expr)->get_symbol()->get_name();
                    CircularArrMap::iterator it2 = transMap.find(rName);
                    if (it2 != transMap.end()) 
                        transMap[lName] = it2->second;
                    else transMap[lName] = rName;
                }
            }
        }
        if (isSgForStatement(*iter)) /* dfs */ 
            getCircularArrays(cirArrMap, isSgForStatement(*iter));
    }

    /* summary */
    FOREACH(CircularArrMap, transMap, iter) {
        CircularArrMap::iterator it = transMap.find(iter->second);
        if (it != transMap.end() && transMap.find(it->second) == iter) {
            cirArrMap[iter->first] = iter->second;
            cirArrMap[iter->second] = iter->first;
        }
    }
}

/* 
 * get useful stencil information for physis
 * precondition: the for nest is a Stencil 
 */
void getStencilInfo(StencilInfo* info, 
        SgNode* node, StructFieldMap &structFieldMap) {
    /* #1 get info->st_node */
    info->st_nodes.push_back(node);

    /* get loop range info list and the innest for loop */
    RangeList loopRangeList; 
    SgNode *outerFor = node, *innerFor;
    getLoopRange(loopRangeList, outerFor);
    while ((innerFor = haveInnerFor(outerFor))) {
        getLoopRange(loopRangeList, innerFor);
        outerFor = innerFor;
    }
   
    /* get loop body of the innest for loop */
    SgNodeList exprList;
    getLoopBody(exprList, outerFor);
    
    DomainInfo domInfo;
    ArrayInfoList arrayInfoList;
    info->st_type = StencilInfo::ST_TYPE_UNSET; 
    FOREACH(SgNodeList, exprList, iter) {
        SgExprStatement* exprStmt = isSgExprStatement(*iter);
        assert(exprStmt);
        SgAssignOp* assignNode = isSgAssignOp(exprStmt->get_expression());
        if (!assignNode) continue; /* TODO: other expressions */

        /* #5 get info->st_type */
        SgType* type = assignNode->get_type();
        if (isSgTypedefType(type))
            type = isSgTypedefType(type)->get_base_type();
        if (isSgTypeFloat(type)) {
            if (info->st_type != StencilInfo::ST_TYPE_UNSET && 
                    info->st_type != StencilInfo::ST_TYPE_FLOAT) {
                info->st_type = StencilInfo::ST_TYPE_UNKNOWN;
            } else info->st_type = StencilInfo::ST_TYPE_FLOAT;
        }
        else if (isSgTypeDouble(type)) {
            if (info->st_type != StencilInfo::ST_TYPE_UNSET && 
                    info->st_type != StencilInfo::ST_TYPE_DOUBLE) {
                info->st_type = StencilInfo::ST_TYPE_UNKNOWN;
            } else info->st_type = StencilInfo::ST_TYPE_DOUBLE;
        } else info->st_type = StencilInfo::ST_TYPE_UNKNOWN;

        SgExpression *lhs = assignNode->get_lhs_operand();
        SgExpression *rhs = assignNode->get_rhs_operand();
        
        ArrayInfoList lhArrayInfoList, rhArrayInfoList;
        getArrayInfo4Expr(lhArrayInfoList, (SgNode* )lhs);
        getArrayInfo4Expr(rhArrayInfoList, (SgNode* )rhs);
        arrayInfoList.insert(arrayInfoList.end(), 
                lhArrayInfoList.begin(), lhArrayInfoList.end());
        arrayInfoList.insert(arrayInfoList.end(), 
                rhArrayInfoList.begin(), rhArrayInfoList.end());

        /* #2 get info->st_domains.d_vars */
        getVarInfo(domInfo.d_vars, (SgNode* )lhs);
        getVarInfo(domInfo.d_vars, (SgNode* )rhs);
    }
    
    /* #3 get info->st_dim */
    info->st_dim = 0;
    FOREACH(ArrayInfoList, arrayInfoList, iter) {
        ArrayInfo arrayInfo = *iter;
        int dim = arrayInfo.second.size();
        if (!info->st_dim)
            info->st_dim = dim;
        else { 
            /* all array in a stencil must have the same dimension */;
            assert(dim == info->st_dim);
        } 
    }
    
    int baseIdx = (int)loopRangeList.size() > info->st_dim;
    
    /* #8 get info->st_dimMaxSize */
    bool isDimMaxSizeSet = false;
    FOREACH(ArrayInfoList, arrayInfoList, iter) {
        SgNode* arrNode = iter->first; 

        /* get var name */
        SgName varName;
        if (isSgVarRefExp(arrNode))
            varName = isSgVarRefExp(arrNode)->get_symbol()->get_name();
        SgArrowExp* arrowExp = isSgArrowExp(arrNode);
        if (arrowExp) {
            varName = getName4ArrowExp(arrowExp);
        }

        SgNodeList shape = getShape(arrNode, structFieldMap);
        
        /* TODO: get shape of the dynamically allocated arrays */
        /*
        if (shape.empty()) {
            for (int i = baseIdx; i < (int)loopRangeList.size(); i++) 
                shape.push_back(loopRangeList[i].second.second);
        }
        */

        /* #9 get info->st_domains.d_arrays */
        domInfo.d_arrays[varName] = shape;

        if (!isDimMaxSizeSet) {
            /* todo get dim MAX size */
            info->st_dimMaxSize.insert(info->st_dimMaxSize.end(), 
                    shape.begin(), shape.end());
            isDimMaxSizeSet = true;
        }
    }
    
    /* #4 get info->st_arrays */
    info->st_arrays = domInfo.d_arrays; 

    /* #6 get info->st_domains */
    domInfo.d_range.assign(loopRangeList.begin() + baseIdx, loopRangeList.end());
    domInfo.d_exprs = exprList;
    info->st_domains.push_back(domInfo);

    /* #7 get info->st_iterate */
    IterateInfo* itInfo = &info->st_iter;
    if (!baseIdx) itInfo->it_times = make_pair(buildIntVal(0), buildIntVal(1));
    else itInfo->it_times = loopRangeList[0].second;
    CircularArrMap cirArrMap;
    getCircularArrays(cirArrMap, isSgForStatement(node));
    itInfo->it_cirarrs = cirArrMap; 
    itInfo->it_circular = !cirArrMap.empty();
}

/* merge stencil info src into dst */
void mergeStencilInfo(StencilInfo &dst, StencilInfo &src) {
    dst.st_domains.insert(dst.st_domains.end(), 
            src.st_domains.begin(), src.st_domains.end());
    
    /* TODO: st_iter, st_dimMaxSize */

    /* st_nodes */
    dst.st_nodes.insert(dst.st_nodes.end(), 
            src.st_nodes.begin(), src.st_nodes.end());

    /* st_arrays */
    FOREACH(ArrayRangeMap, src.st_arrays, iter) {
        dst.st_arrays[iter->first] = iter->second;
    }

    /* st_iter.it_cirarrs */
    dst.st_iter.it_cirarrs.insert(
            src.st_iter.it_cirarrs.begin(), src.st_iter.it_cirarrs.end());
}

/* merge stencils */
void mergeStencilInfoList(StencilInfoList &stencilInfoList) {
    ForStencilInfoMap forStencilInfoMap;
    FOREACH(StencilInfoList, stencilInfoList, iter) {
        forStencilInfoMap[iter->st_nodes[0]] = *iter; 
    }

    /* 
     * if the for nest after current for nest is also a stencil,
     * and the st_dim, st_type are the same
     * TODO: the it_times are the same
     */
    FOREACH(ForStencilInfoMap, forStencilInfoMap, iter) {
        SgStatement* curForNest = isSgStatement(iter->first);
        while (true) {
            SgStatement* nextForNest = getNextStatement(curForNest);
            ForStencilInfoMap::iterator it = forStencilInfoMap.find(nextForNest);
            if (it != forStencilInfoMap.end() && 
                    it->second.st_dim == iter->second.st_dim && 
                    it->second.st_type == iter->second.st_type) {
                mergeStencilInfo(iter->second, it->second); 
                curForNest = nextForNest;
                forStencilInfoMap.erase(it);
            } else break;
        }
    }

    stencilInfoList.clear();
    FOREACH(ForStencilInfoMap, forStencilInfoMap, iter) {
        stencilInfoList.push_back(iter->second); 
    }
}


/* judge wether the call-site is in a initializer of a struct */
SgType* isInitOfStruct(SgNode* allocSite) {
    SgNode* parent = allocSite->get_parent();
    if (isSgCastExp(parent)) 
        parent = parent->get_parent();
    if (!isSgAssignInitializer(parent)) return NULL;
    
    /* xxx = calloc() */
    parent = parent->get_parent();
    SgDesignatedInitializer* designatedIniter 
        = isSgDesignatedInitializer(parent);
    if (!designatedIniter) return NULL;

    /* .xxx = calloc() */
    parent = parent->get_parent();
    if (!isSgExprListExp(parent)) return NULL;
    
    /* = {..., .xxx = calloc(), ...}*/
    parent = parent->get_parent();
    if (!isSgAggregateInitializer(parent)) return NULL;

    /* yyy = {..., .xxx = calloc(), ...} */
    parent = parent->get_parent();
    SgInitializedName *initName = isSgInitializedName(parent);
    if (!initName) return NULL;
    
    /* struct ZZZ yyy = {..., .xxx = calloc(), ...} */
    SgType* type = initName->get_type();
    if (!isStructType(type)) return NULL;

    return type; 
}

/* .field = calloc(), get field name */
SgName getInitField(SgNode* allocSite) {
    SgName null = SgName();

    /* find designated initializer */
    SgNode* parent = allocSite->get_parent();
    while (!isSgDesignatedInitializer(parent))
        parent = parent->get_parent();

    SgDesignatedInitializer* designatedIniter 
        = isSgDesignatedInitializer(parent);
    SgExpressionPtrList expList 
        = designatedIniter->get_designatorList()->get_expressions();
    if (!isSgVarRefExp(expList[0])) return null;
    
    return isSgVarRefExp(expList[0])->get_symbol()->get_name();
}

/* .field = symbol, get field name */
SgName matchField(SgVariableSymbol* symbol, SgNode* allocSite) {
    
    /* find exprListExp of AggregateInitializer */
    SgNode* parent = allocSite->get_parent();
    while (!isSgExprListExp(parent))
        parent = parent->get_parent();
    SgExpressionPtrList fieldList = isSgExprListExp(parent)->get_expressions();
    FOREACH(SgExpressionPtrList, fieldList, iter) {
        SgDesignatedInitializer* designatedIniter 
            = isSgDesignatedInitializer(*iter);
        SgAssignInitializer* assignIniter 
            = isSgAssignInitializer(designatedIniter->get_memberInit());
        if (!assignIniter) continue;

        /* .xxx = something simple */
        SgExpression* rhs = assignIniter->get_operand();
        if (!isSgVarRefExp(rhs)) continue;
        
        /* .xxx = yyy */
        if (isSgVarRefExp(rhs)->get_symbol() == symbol) {
            SgExpressionPtrList expList 
                = designatedIniter->get_designatorList()->get_expressions();
            return isSgVarRefExp(expList[0])->get_symbol()->get_name();
        }
    }
    return SgName();
} 

void getFieldList(FieldList &fieldList, SgNode* allocSite) {
    SgFunctionCallExp *funcCallExp = isSgFunctionCallExp(allocSite);
    SgExpressionPtrList argList = funcCallExp->get_args()->get_expressions();
    
    /* 1D: return argument 1 */
    if (isSgVarRefExp(argList[0])) {
        SgVariableSymbol* symbol = isSgVarRefExp(argList[0])->get_symbol();        
        SgName field = matchField(symbol, allocSite);
        fieldList.push_back(field);
    } else {
        /* other dimessions: return sizeof([][]) */
        SgType* opType = isSgSizeOfOp(argList[1])->get_operand_type();
        SgExpressionList expList 
            = get_C_array_dimensions(*isSgArrayType(opType));
        FOREACH(SgExpressionList, expList, iter) {
            SgVariableSymbol* symbol = isSgVarRefExp(*iter)->get_symbol();
            SgName field = matchField(symbol, allocSite);
            fieldList.push_back(field);
        }
    }
}

void fillStructFieldLength(SgProject* project, StructFieldMap &structFieldMap) {
    /* find all call-site to "alloc" functions */
    SgNodeSet allocSites;
    getAllocSite(allocSites, project);
   
    /* for each alloc-site, fill the structure */
    FOREACH(SgNodeSet, allocSites, iter) {
        /* must be initialization of a struct */
        SgType* structType = isInitOfStruct(*iter);
        if (!structType) continue;
        
        /* get field name */
        SgName field = getInitField(*iter);
        if (field.is_null()) continue;
       
        FieldList fieldList;
        getFieldList(fieldList, *iter);
        
        if (structFieldMap.find(structType) == structFieldMap.end()) {
            FieldRangeMap fieldRangeMap;
            fieldRangeMap[field] = fieldList;
            structFieldMap[structType] = fieldRangeMap;
        } else {
            structFieldMap[structType][field] = fieldList; 
        }
    }
    
    dbg << GAP << "struct field map\n" << GAP;
    FOREACH(StructFieldMap, structFieldMap, iter) {
        dbg << iter->first->unparseToString() << "\n";
        FOREACH(FieldRangeMap, iter->second, it) {
            dbg << "\t" << it->first.getString() << ": ";
            FOREACH(FieldList, it->second, it2) 
                dbg << (*it2).getString() << ", "; 
            dbg << "\n";
        }
        dbg << "\n";
    }
}

/* get stencil info for stencil for nest */
void getStencilInfoList(StencilInfoList &stencilInfoList, 
        SgNodeList &stencilNodeList, StructFieldMap &structFieldMap) {

    /* collect stencil info */ 
    FOREACH(SgNodeList, stencilNodeList, iter) {
        StencilInfo info;
        getStencilInfo(&info, *iter, structFieldMap);
        stencilInfoList.push_back(info);
    }
    
    dbg << GAP << "stencil info BEFORE MERGE\n" << GAP;
    int stencilId = 0;
    FOREACH(StencilInfoList, stencilInfoList, iter) {    
        printStencilInfo(&(*iter), ++stencilId);
        dbg << "\n";
    }
    dbg << "# TOTALLY " << stencilInfoList.size() << " STENCILS\n\n";
   
    mergeStencilInfoList(stencilInfoList);
    
    dbg << GAP << "stencil info AFTER MERGE\n" << GAP;
    stencilId = 0;
    FOREACH(StencilInfoList, stencilInfoList, iter) {    
        printStencilInfo(&(*iter), ++stencilId);
        dbg << "\n";
    }
    dbg << "# TOTALLY " << stencilInfoList.size() << " STENCILS\n\n";
}
	
// ##############################
//          modify code
// ##############################
/* add #include "physis/physis.h" */
/*
void myInsertHeader(SgProject* project) {
	SgGlobal* globalscope = getFirstGlobalScope(project);
	insertHeader("physis/physis.h", PreprocessingInfo::before, 
					false, globalscope);
}
*/

SgFunctionDeclaration* getFuncDecl4Node(SgNode *node) {
    SgNode* scope = node;
    while (scope) {
        scope = scope->get_parent();
        SgFunctionDeclaration* funcDecl = isSgFunctionDeclaration(scope);
        if (funcDecl) return funcDecl;
    }
    return NULL;
}

/* 
 * check wether the scope of stencil 
 * have variable int argc and char* argv[],
 * if not, pass them as global variable
 */
string dealWithArgc(string argcName, 
        SgStatement* stmt, SgScopeStatement* globalScope) {
    
    dbg << "deal with argc\n";

    SgVariableSymbol* symbolArgc 
        = lookupVariableSymbolInParentScopes(SgName("argc"), stmt->get_scope());
    /* if not, append them to the parameter list */
    if (!symbolArgc) { 
        argcName = "globalArgc";
        /* 
         * add "extern int globalArgc;" 
         * to current file if globalArgc not exist 
         */
        symbolArgc = lookupVariableSymbolInParentScopes(
                SgName(argcName), stmt->get_scope());
        if (!symbolArgc) {
            /* insert "extern int globalArgc;" */
	        SgVariableDeclaration* varDecl 
	            = buildVariableDeclaration(
	                    argcName, 
	                    buildIntType(), 
	                    NULL,
	                    globalScope);
	        setExtern(varDecl);
	        insertStatement(getFirstStatement(globalScope), varDecl);
        }
    }
    return argcName;
}

/* 
 * check wether the scope of stencil 
 * have variable int argc and char* argv[],
 * if not, pass them as global variable
 */
string dealWithArgv(string argvName, 
        SgStatement* stmt, SgScopeStatement* globalScope) {
    
    dbg << "deal with argv\n";

    SgVariableSymbol* symbolArgv 
        = lookupVariableSymbolInParentScopes(SgName("argv"), stmt->get_scope());
    /* if not, append them to the parameter list */
    if (!symbolArgv) { 
        argvName = "globalArgv";
        /*
         * add "extern int globalArgv;" 
         * to current file if globalArgv not exist
         */
        symbolArgv = lookupVariableSymbolInParentScopes(
                SgName(argvName), stmt->get_scope());
        if (!symbolArgv) {
            /* insert "extern char** globalArgv;" */
	        SgVariableDeclaration* varDecl 
	            = buildVariableDeclaration(
                        argvName,
	                    buildPointerType(buildPointerType(buildCharType())), 
	                    NULL,
	                    globalScope);
	        setExtern(varDecl);
	        insertStatement(getFirstStatement(globalScope), varDecl);
        }
    }
    return argvName;
}

/* insert void PSInit(int *argc, char ***argv, int grid_num_dims, ...) */
SgStatement* buildPSInit(StencilInfo *info, 
        SgStatement* stmt, SgScopeStatement* globalScope) {
    string argcName = dealWithArgc("argc", stmt, globalScope);
    
    string argvName = dealWithArgv("argv", stmt, globalScope);

    dbg << "insert PSInit()\n\n";
    
    SgFunctionDeclaration *funcDecl = getFuncDecl4Node(info->st_nodes[0]);
    // SgFunctionParameterList* paramList = funcDecl->get_parameterList();
    SgFunctionDefinition* funcDef = funcDecl->get_definition();
    
    /* build call statement */
    SgType* returnType = buildVoidType();
    SgExprListExp* argList = buildExprListExp();
    SgAddressOfOp* argArgc = buildAddressOfOp(buildVarRefExp(argcName, funcDef));
    SgAddressOfOp* argArgv = buildAddressOfOp(buildVarRefExp(argvName, funcDef));
    SgIntVal* argDim = buildIntVal(info->st_dim);
    appendExpression(argList, argArgc);
    appendExpression(argList, argArgv);
    appendExpression(argList, argDim);
    FOREACH(SgNodeList, info->st_dimMaxSize, iter) {
        appendExpression(argList, copyExpression(isSgExpression(*iter)));
    }
    SgExprStatement* callPSInit = buildFunctionCallStmt(SgName("PSInit"), 
            returnType, argList, stmt->get_scope());
    
    /* insert call statement */
    insertStatement(stmt, callPSInit);

    return callPSInit;
}

/* insert void PSFinalize() */
SgStatement* buildPSFinalize(SgStatement *stmt) {
    
    dbg << "insert PSFinalize()\n\n"; 
    
    /* build call statement */
    SgType* returnType = buildVoidType();
    SgExprStatement* callPSFinalize = buildFunctionCallStmt(SgName("PSFinalize"), 
            returnType, NULL, stmt->get_scope());
    
    /* insert call statement */
    insertStatement(stmt, callPSFinalize, false);
    
    return callPSFinalize;
}

/* fill <dim> and <type> for PSGrid<dim><type>New */
string getGridTypeName(int dim, StencilInfo::ST_TYPE type) {
    string result = "PSGrid";
    result += boost::lexical_cast<string>(dim) + "D";
    
    switch (type) {
        case StencilInfo::ST_TYPE_DOUBLE:
            result += "Double";
            break;
        case StencilInfo::ST_TYPE_FLOAT:
            result += "Float";
            break;
        default:
            dbg << "Something wrong with stencil info\n";
            break;
    }
    return result;
}

/* 
 * generate variable name in kernel, 
 * e.g. fdtd->dx to fdtd_dx, dx to dx 
 */
string generateParamName(string varName) {
    string::size_type pos = varName.rfind("->");
    if (pos != varName.npos)
        varName = varName.substr(0, pos) 
            + varName.substr(pos + 2, varName.length() - pos - 2);
    return varName;
}

/* 
 * generate gird name for array, 
 * e.g. fdtd->ez to gez, arr to garr
 */
string generateGridName(string arrName) {
    string::size_type pos = arrName.rfind("->");
    if (pos != arrName.npos)
        arrName = arrName.substr(pos + 2, arrName.length() - pos - 2);
    return "g" + arrName;
}

/* insert PSGrid<dim><type>New(); */
SgStatement* buildPSGridNew(SgType* gridType, string funcName, 
        string arrName, SgNodeList &dimInfo, SgStatement* stmt) {
    
    dbg << "insert PSGrid<dim><type>New() for " << arrName << "\n"; 
    
    /* build call statement */
    SgExprListExp* argList = buildExprListExp();
    FOREACH(SgNodeList, dimInfo, iter) {
        appendExpression(argList, copyExpression(isSgExpression(*iter)));
    }
    SgFunctionCallExp* callPSGridNew 
        = buildFunctionCallExp(SgName(funcName), 
            gridType, argList, stmt->get_scope());
    /* build variable declaration */
    SgAssignInitializer* initer
        = buildAssignInitializer(callPSGridNew, gridType); 
    string gridName = generateGridName(arrName);
    SgVariableDeclaration* varDecl 
        = buildVariableDeclaration(SgName(gridName), gridType, 
                initer, stmt->get_scope());

    /* insert call statement */
    insertStatement(stmt, varDecl, false);
    
    return varDecl;
}

/* insert void PSGridCopyin(void* dst, const void* src); */
SgStatement* buildPSGridCopyin(string arrName, SgStatement* stmt) {
    
    dbg << "insert PSGridCopyin() for " << arrName << "\n"; 
    
    /* build call statement */
    SgType* returnType = buildVoidType();
    SgVarRefExp* argGArr 
        = buildVarRefExp(generateGridName(arrName), stmt->get_scope());
    SgVarRefExp* argArr = buildVarRefExp(arrName, stmt->get_scope());
    SgExprListExp* argList = buildExprListExp();
    appendExpression(argList, argGArr);
    appendExpression(argList, argArr);
    SgExprStatement* callPSGridCopyin 
        = buildFunctionCallStmt(SgName("PSGridCopyin"), 
            returnType, argList, stmt->get_scope());

    /* insert call statement */
    insertStatement(stmt, callPSGridCopyin, false);
    
    return callPSGridCopyin;
}

/* insert void PSGridFree(void* p); */
SgStatement* buildPSGridFree(string arrName, SgStatement* stmt) {
    
    dbg << "insert PSGridFree() for " << arrName << "\n"; 
    
    /* build call statement */
    SgType* returnType = buildVoidType();
    SgVarRefExp* argGArr 
        = buildVarRefExp(generateGridName(arrName), stmt->get_scope());
    SgExprListExp* argList = buildExprListExp();
    appendExpression(argList, argGArr);
    SgExprStatement* callPSGridFree 
        = buildFunctionCallStmt(SgName("PSGridFree"), 
            returnType, argList, stmt->get_scope());


    /* insert call statement */
    insertStatement(stmt, callPSGridFree, false);
    
    return callPSGridFree;
}

/* fill <dim> for PSDomain<dim>New */
string getDomainTypeName(int dim) {
    string result = "PSDomain";
    result += boost::lexical_cast<string>(dim) + "D";
    return result;
}

/* insert PSDomain<dim>New(); */
SgStatement* buildPSDomainNew(SgType* domainType, string funcName, 
        RangeList &rangeList, SgName domainName, SgStatement* stmt) {

    dbg << "insert PSDomain<dim>New() for " << domainName.getString() << " \n"; 
    
    /* build call statement */
    SgExprListExp* argList = buildExprListExp();
    FOREACH(RangeList, rangeList, iter) {
        myRange range = (*iter).second;
        appendExpression(argList, copyExpression(isSgExpression(range.first)));
        appendExpression(argList, copyExpression(isSgExpression(range.second)));
    }
    SgFunctionCallExp* callPSDomainNew 
        = buildFunctionCallExp(SgName(funcName), 
            domainType, argList, stmt->get_scope());
    
    /* build variable declaration */
    SgAssignInitializer* initer
        = buildAssignInitializer(callPSDomainNew, domainType); 
    SgVariableDeclaration* varDecl 
        = buildVariableDeclaration(domainName, domainType, 
                initer, stmt->get_scope());

    /* insert call statement */
    insertStatement(stmt, varDecl, false);
    
    return varDecl;
}

/* build PSStencil PSStencilMap(void *p, ...) */
SgExpression* buildPSStencilMap(DomainInfo &domInfo, 
        SgTypedefTypeMap &typedefMap,  
        string kernelName, string domainName, SgStatement* stmt,
        CircularArrMap* cirArrMap = NULL, bool swap = false) {
    /* return type */
    SgType* returnType = typedefMap[SgName("PSStencil")];  
    
    /* argument list */
    SgExprListExp* argList = buildExprListExp();
    SgVarRefExp* kernel = buildVarRefExp(kernelName, stmt->get_scope());
    SgVarRefExp* domain = buildVarRefExp(domainName, stmt->get_scope());
    appendExpression(argList, kernel);
    appendExpression(argList, domain);
    /* pass arrays as arguments */
    FOREACH(ArrayRangeMap, domInfo.d_arrays, iter) {
        string arrName = iter->first.getString();
        string gridName = generateGridName(arrName);
        if (swap) {
            CircularArrMap::iterator it = cirArrMap->find(iter->first);
            if (it != cirArrMap->end()) {
                gridName = generateGridName(it->second.getString());
            }
        } 
        SgVarRefExp* arrRef = buildVarRefExp(gridName, stmt->get_scope());
        appendExpression(argList, arrRef);
    }
    /* pass vars as arguments */
    FOREACH(VarInfoMap, domInfo.d_vars, iter) {
        SgVarRefExp* varRef = buildVarRefExp(iter->first, stmt->get_scope());
        appendExpression(argList, varRef);
    } 
    
    /* call expression */
    return buildFunctionCallExp(SgName("PSStencilMap"), 
            returnType, argList, stmt->get_scope());
}


string generateDomainName(int stencilId, int kernelId) {
    string suffix = "_S" + boost::lexical_cast<string>(stencilId) 
        + "D" + boost::lexical_cast<string>(kernelId); 
    return "domain" + suffix;
}

string generateKernelName(int stencilId, int kernelId) {
    string suffix = "_S" + boost::lexical_cast<string>(stencilId) 
        + "K" + boost::lexical_cast<string>(kernelId); 
    return "kernel" + suffix;
}

/* insert void PSStencilRun(PSStencil, ...) */
SgStatement* buildPSStencilRun(StencilInfo* info, SgTypedefTypeMap &typedefMap, 
        int stencilId, CircularArrMap &cirArrMap, SgStatement* stmt) {
    
    dbg << "insert PSStencilRun()\n\n"; 
    
    /* return type */   
    SgType* returnType = buildVoidType();

    /* argument list */
    SgExprListExp* argList = buildExprListExp();
    if (!info->st_iter.it_circular) {
        /* pass the return value of PSStencilMap() as argument */
        int kernelId = 0;
        FOREACH(DomainInfoList, info->st_domains, iter) {
            kernelId++;
            string kernelName = generateKernelName(stencilId, kernelId);
            string domainName = generateDomainName(stencilId, kernelId);
            SgExpression* exp 
                = buildPSStencilMap(*iter, typedefMap, kernelName, domainName, stmt);
            appendExpression(argList, exp);
        }
        /* pass iterate times as argument */
        SgExpression *lRange = isSgExpression(info->st_iter.it_times.first);
        SgExpression *rRange = isSgExpression(info->st_iter.it_times.second);
        appendExpression(argList, 
                buildSubtractOp(copyExpression(rRange), 
                    copyExpression(lRange)));
    } else {
        /* pass the return value of PSStencilMap() as argument */
        int kernelId = 0;
        FOREACH(DomainInfoList, info->st_domains, iter) {
            kernelId++;
            string kernelName = generateKernelName(stencilId, kernelId);
            string domainName = generateDomainName(stencilId, kernelId);
            SgExpression* exp1
                = buildPSStencilMap(*iter, typedefMap, kernelName, domainName, stmt, 
                        &cirArrMap);
            SgExpression* exp2
                = buildPSStencilMap(*iter, typedefMap, kernelName, domainName, stmt, 
                        &cirArrMap, true);
            appendExpression(argList, exp1);
            appendExpression(argList, exp2);
        }
        /* pass iterate times as argument */
        SgExpression *lRange = isSgExpression(info->st_iter.it_times.first);
        SgExpression *rRange = isSgExpression(info->st_iter.it_times.second);
        /* (rRange - lRange + 1) / 2 */
        SgSubtractOp *subOp 
            = buildSubtractOp(copyExpression(rRange), copyExpression(lRange));
        SgAddOp *addOp = buildAddOp(subOp, buildIntVal(1));
        SgDivideOp *divOp = buildDivideOp(addOp, buildIntVal(2));
        appendExpression(argList, divOp);
    }
   
    /* call statement */
    SgExprStatement* callPSStencilRun 
        = buildFunctionCallStmt(SgName("PSStencilRun"), 
            returnType, argList, stmt->get_scope());
   
    /* insert call statement */
    insertStatement(stmt, callPSStencilRun, false);

    return callPSStencilRun;
}

/* insert void PSGridCopyout(void* p, const void* dst_arr); */
SgStatement* buildPSGridCopyout(myRange &it_times, string arrName, 
        CircularArrMap& cirArrMap, SgStatement* stmt) {
    
    dbg << "insert PSGridCopyout() for " << arrName << "\n"; 
    
    /* not using circular array */
    if (cirArrMap.find(SgName(arrName)) == cirArrMap.end()) {
	    /* build call statement */
	    SgType* returnType = buildVoidType();
	    SgVarRefExp* argGrid 
            = buildVarRefExp(generateGridName(arrName), stmt->get_scope());
	    SgVarRefExp* argArr = buildVarRefExp(arrName, stmt->get_scope());
	    SgExprListExp* argList = buildExprListExp();
	    appendExpression(argList, argGrid);
	    appendExpression(argList, argArr);
	    SgExprStatement* callPSGridCopyout 
	        = buildFunctionCallStmt(SgName("PSGridCopyout"), 
	            returnType, argList, stmt->get_scope());
	
	    /* insert call statement */
	    insertStatement(stmt, callPSGridCopyout, false);
	    
	    return callPSGridCopyout;
    } /* else: using circular array */

    /* build if statement */
    
    /* conditional */
    SgModOp* conditional = buildModOp(
            buildSubtractOp(
                copyExpression(isSgExpression(it_times.second)), 
                copyExpression(isSgExpression(it_times.first))),
            buildIntVal(2));
    
    /* true body */
    SgType* returnType = buildVoidType();
    /* move u <-- guNext or uNext <-- gu */
    SgVarRefExp* argGrid 
        = buildVarRefExp(
                generateGridName(cirArrMap[SgName(arrName)].getString()), 
                stmt->get_scope());
    SgVarRefExp* argArr = buildVarRefExp(arrName, stmt->get_scope());
    SgExprListExp* argList = buildExprListExp();
    appendExpression(argList, argGrid);
    appendExpression(argList, argArr);
    SgExprStatement* trueBody 
        = buildFunctionCallStmt(SgName("PSGridCopyout"), 
            returnType, argList, stmt->get_scope());
    
    /* false body */
    /* move u <-- gu or uNext <-- guNext */
    argGrid = buildVarRefExp(generateGridName(arrName), stmt->get_scope());
    argArr = buildVarRefExp(arrName, stmt->get_scope());
    argList = buildExprListExp();
    appendExpression(argList, argGrid);
    appendExpression(argList, argArr);
    SgExprStatement* falseBody 
        = buildFunctionCallStmt(SgName("PSGridCopyout"), 
            returnType, argList, stmt->get_scope());

    SgIfStmt *ifStmt = buildIfStmt(conditional, trueBody, falseBody);
    
    /* insert if statement */
    insertStatement(stmt, ifStmt, false);
    
    return ifStmt;
}

/* my build function parameter list */
SgFunctionParameterList* buildFuncParList(int dim, DomainInfo &domInfo, 
        StrList &idList, SgType* gridType) {
    SgFunctionParameterList* parList = buildFunctionParameterList();
    
    /* indexs */
    for (int i = 0; i < dim; i++)
        appendArg(parList, 
                buildInitializedName(
                    idList[i], 
                    buildConstType(buildIntType()))
                );

    /* girds */
    FOREACH(ArrayRangeMap, domInfo.d_arrays, iter) {
        appendArg(parList, 
                buildInitializedName(
                    generateGridName(iter->first.getString()), 
                    gridType)
                );    
    }
    /* variables */
    FOREACH(VarInfoMap, domInfo.d_vars, iter) {
        appendArg(parList, buildInitializedName(
                    generateParamName(iter->first.getString()), 
                    iter->second)
                ); 
    }
    return parList;
}

/* copy the "i +- c" expression and left out all type casts */
SgExpression* rebuildIndex(SgNode* node, SgScopeStatement* scope) {
    SgExpression* result = NULL;
    SgExpression* exp = isSgExpression(node);
    if (isSgAddOp(exp)) {
        SgExpression* lhs = isSgAddOp(exp)->get_lhs_operand();
        SgExpression* rhs = isSgAddOp(exp)->get_rhs_operand();
        SgName varName;
        if (isSgVarRefExp(lhs)) {
            varName = isSgVarRefExp(lhs)->get_symbol()->get_name();
        }
        int intVal = 0;
        if (isSgIntVal(rhs)) {
            intVal = isSgIntVal(rhs)->get_value();
        } else if (isSgCastExp(rhs)) {
            intVal = isSgIntVal(isSgCastExp(rhs)->get_operand())->get_value();
        }
        result = buildAddOp(
                buildVarRefExp(varName, scope), 
                buildIntVal(intVal));
    }
    else if (isSgSubtractOp(exp)) {
        SgExpression* lhs = isSgSubtractOp(exp)->get_lhs_operand();
        SgExpression* rhs = isSgSubtractOp(exp)->get_rhs_operand();
        SgName varName;
        if (isSgVarRefExp(lhs)) {
            varName = isSgVarRefExp(lhs)->get_symbol()->get_name();
        }
        int intVal = 0;
        if (isSgIntVal(rhs)) {
            intVal = isSgIntVal(rhs)->get_value();
        } else if (isSgCastExp(rhs)) {
            intVal = isSgIntVal(isSgCastExp(rhs)->get_operand())->get_value();
        }
        result = buildSubtractOp(
                buildVarRefExp(varName, scope),
                buildIntVal(intVal)
                );
    
    } else if (isSgVarRefExp(exp)) {
        SgName varName = isSgVarRefExp(exp)->get_symbol()->get_name();
        result = buildVarRefExp(varName, scope);
    }
    return result;
}

/* change the array reference in stencil statement into PSGridGet() */
void changeArr2Grid(SgExpression* node, SgFunctionDefinition* funcDef) { 
    SgPntrArrRefExp* arrRef = isSgPntrArrRefExp(node);
    if (arrRef) {
        /* get array info */
        ArrayInfo info;
        getArrayInfo4Arr(info, arrRef);
        
        /* build PSGridGet() */
        SgType* returnType = arrRef->get_type();
        SgExprListExp* argList = buildExprListExp();
        string arrName;
        if (isSgVarRefExp(info.first))
            arrName = isSgVarRefExp(info.first)->get_symbol()
                ->get_name().getString();
        else if (isSgArrowExp(info.first))
            arrName = getName4ArrowExp(isSgArrowExp(info.first));
        appendExpression(argList, 
                buildVarRefExp(
                    generateGridName(arrName),
                    funcDef)
                );
        for (int i = info.second.size() - 1; i >=0 ; i--) {
            SgExpression* exp = rebuildIndex(info.second[i], funcDef);
            appendExpression(argList, exp);
            /* replace variable reference */
            /*
            SgNodeList varRefList = NodeQuery::querySubTree(info.second[i], V_SgVarRefExp);
            FOREACH(SgNodeList, varRefList, iter) {
                SgVarRefExp* varRef = isSgVarRefExp(*iter);
                SgVariableSymbol* symbol 
                    = lookupVariableSymbolInParentScopes(varRef->get_symbol()->get_name(), funcDef);
                varRef->set_symbol(symbol); 
            }
            appendExpression(argList, 
                    copyExpression(isSgExpression(info.second[i])));
            */
        }

        SgExpression* callPSGridGet 
            = buildFunctionCallExp(
                    SgName("PSGridGet"), 
                    returnType, 
                    argList, 
                    funcDef);

        replaceExpression(arrRef, callPSGridGet);
        return;
    }

    SgNodeList children = node->get_traversalSuccessorContainer();
    FOREACH(SgNodeList, children, iter) {
        changeArr2Grid(isSgExpression(*iter), funcDef);
    }
}

/* change the variable reference in stencil statement into passed parameters */
void changeVar(SgExpression* node, SgFunctionDefinition* funcDef) { 
    if (isSgPntrArrRefExp(node)) return;

    if (isSgArrowExp(node)) {
        string varName = getName4ArrowExp(isSgArrowExp(node));
        string paramName = generateParamName(varName);
        replaceExpression(node, buildVarRefExp(paramName, funcDef));
        return;
    }

    if (isSgVarRefExp(node)) {
        string varName = isSgVarRefExp(node)->get_symbol()
            ->get_name().getString();
        string paramName = generateParamName(varName);
        replaceExpression(node, buildVarRefExp(paramName, funcDef));
        return;
    }
    
    SgNodeList children = node->get_traversalSuccessorContainer();
    FOREACH(SgNodeList, children, iter) {
        changeVar(isSgExpression(*iter), funcDef);
    }
}

/* insert PSGridEmit() */
SgStatement* buildPSGridEmit(string arrName, SgExpression* rhs, 
        SgFunctionDefinition* funcDef) {
    
    dbg << "insert PSGridEmit() for " << arrName << "\n"; 
    
    /* return type */
    SgType* returnType = buildVoidType();
    
    /* argument list */
    SgExprListExp* argList = buildExprListExp();
    appendExpression(argList, 
            buildVarRefExp(
                generateGridName(arrName), 
                funcDef)
            );
    appendExpression(argList, rhs);

    return buildFunctionCallStmt(
            SgName("PSGridEmit"), 
            returnType, 
            argList, 
            funcDef);
}



/* modify each statement in stencil for nest */
SgStatement* modifyStatement(SgNode* node, SgFunctionDefinition* funcDef) {
    SgExprStatement* exprStmt = isSgExprStatement(node);
    SgAssignOp *assignOp = isSgAssignOp(exprStmt->get_expression());
    assert(assignOp);

    /* modify right hand side */
    SgExpression* rhs = assignOp->get_rhs_operand();
    changeArr2Grid(rhs, funcDef);
    changeVar(rhs, funcDef);
    
    /* modify left hand side */
    SgExpression* lhs = assignOp->get_lhs_operand();
    ArrayInfo info;
    getArrayInfo4Arr(info, lhs);
   
    string arrName;
    if (isSgVarRefExp(info.first)) 
        arrName = isSgVarRefExp(info.first)->get_symbol()
            ->get_name().getString();
    else if (isSgArrowExp(info.first)) 
        arrName = getName4ArrowExp(isSgArrowExp(info.first));
    
    return copyStatement(
            buildPSGridEmit(arrName, rhs, funcDef));
}

/* insert function body */
void buildFuncBody(SgFunctionDeclaration* funcDecl, DomainInfo &domInfo) {
    SgFunctionDefinition* funcDef = funcDecl->get_definition();
    SgBasicBlock* bb = funcDef->get_body();
   
    FOREACH(SgNodeList, domInfo.d_exprs, iter) 
        appendStatement(modifyStatement(*iter, funcDef), bb);  
}

/* insert kernel function */
void buildKernel(int dim, DomainInfo &domInfo, SgStatement* insertBeforeStmt, 
        SgName kernelName, StrList &idList, 
        SgType* gridType, SgScopeStatement* globalScope) {
    
    dbg << "insert kernel function for " << kernelName.getString() << "\n"; 
    
    /* build function parameter list */
    SgFunctionParameterList* parList 
        = buildFuncParList(dim, domInfo, idList, gridType);
   
    /* build function body */
    SgFunctionDeclaration* funcDecl 
        = buildDefiningFunctionDeclaration(
                kernelName,
                buildVoidType(),
                parList,
                globalScope);
    
    buildFuncBody(funcDecl, domInfo); // TODO: the second parameter should be domain 
    
    /* insert function */
    insertStatement(insertBeforeStmt, funcDecl);
}

/* use physis to represent a single stencil */
void modifyStencil(StencilInfo *info, SgTypedefTypeMap &typedefMap,
        int stencilId, StrList &idList, SgScopeStatement *globalScope) {
    SgStatement* stmt = isSgStatement(info->st_nodes[0]);
    /* build PSInit() */
    stmt = buildPSInit(info, stmt, globalScope);
    
    /* build PSGrid<dim><type>New() */
    string gridTypeName = getGridTypeName(info->st_dim, info->st_type);
    SgType* gridType = typedefMap[SgName(gridTypeName)];
    FOREACH(ArrayRangeMap, info->st_arrays, iter) {
        stmt = buildPSGridNew(gridType, gridTypeName + "New", 
                iter->first.getString(), iter->second, stmt);
    }
    dbg << "\n";
   
    /* build PSGridCopyin() */
    FOREACH(ArrayRangeMap, info->st_arrays, iter) {
        stmt = buildPSGridCopyin(iter->first.getString(), stmt);
    }
    dbg << "\n";
   
    /* build PSDomain<dim>New() */
    /* a stencil may have one or more domains, 
     * each domain is related to a different kernel
     * so the name of domain have format "domain<stencilId>_<kernelId>",
     * and the name of kernel have format "kernel<stencilId>_<kernelId>"
     * */
    string domainTypeName = getDomainTypeName(info->st_dim);
    SgType* domainType = typedefMap[SgName(domainTypeName)];
    int kernelId = 0;
    FOREACH(DomainInfoList, info->st_domains, iter) {
        string domainName = generateDomainName(stencilId, ++kernelId); 
        stmt = buildPSDomainNew(domainType, domainTypeName + "New", 
                iter->d_range, SgName(domainName), stmt);
    }
    dbg << "\n";
   
  /* build kernel function */
    kernelId = 0;
    FOREACH(DomainInfoList, info->st_domains, iter) {
        string kernelName = generateKernelName(stencilId, ++kernelId);
        buildKernel(info->st_dim, *iter, getFuncDecl4Node(info->st_nodes[0]), 
                SgName(kernelName), idList, gridType, globalScope);
        dbg << "\n";
    }
    
    /*
    // TODO: too naive
    CircularArrMap cirArrMap;
    FOREACH(SgNodeList, info->st_nodes, iter) {
        getCircularArrays(cirArrMap, isSgForStatement(*iter));
    }
    */
    
    /* build PSStencilRun */
    stmt = buildPSStencilRun(info, typedefMap, stencilId, 
            info->st_iter.it_cirarrs, stmt);
    
    /* build PSGridCopyout() */
    FOREACH(ArrayRangeMap, info->st_arrays, iter)
        stmt = buildPSGridCopyout(info->st_iter.it_times, iter->first, 
                info->st_iter.it_cirarrs, stmt);
    dbg << "\n";

    /* build PSGridFree() */
    FOREACH(ArrayRangeMap, info->st_arrays, iter) {
        stmt = buildPSGridFree(iter->first, stmt);
    }
    dbg << "\n";

    /* build PSFinalize() */
    stmt = buildPSFinalize(stmt);
    
    FOREACH(SgNodeList, info->st_nodes, iter) {
        removeStatement(isSgStatement(*iter));
        deleteAST(*iter);
    }
} 

/* 
 * store all class type in classTypeMap,
 * and search it with SgName
 */
void buildSgTypedefTypeMap(SgTypedefTypeMap &typedefMap, SgProject* project) {
    SgNodeList result = NodeQuery::querySubTree(project, V_SgTypedefType);
    FOREACH(SgNodeList, result, iter) {
        SgTypedefType *typedefType = isSgTypedefType(*iter);
        typedefMap[typedefType->get_name()] = typedefType;
    }
}

void dealWithArgcArgvInMain(string argcName, string argvName, 
        SgScopeStatement* globalScope) {
    SgFunctionDeclaration* mainDecl = findMain(globalScope);
    if (!mainDecl) return;
    
    /* 
     * insert "int globalArgc;" in globalScope,
     * and "globalArgc = argc;" in main function 
     */
    SgFunctionDefinition* mainDef = mainDecl->get_definition();
    SgVariableDeclaration* varDecl 
        = buildVariableDeclaration(
                argcName, 
                buildIntType(), 
                NULL,
                globalScope);
    insertStatement(getFirstStatement(globalScope), varDecl);
    SgAssignOp* assignOp = buildAssignOp(
            buildVarRefExp(argcName, globalScope),
            buildVarRefExp("argc", mainDef)); 
    mainDef->prepend_statement(buildExprStatement(assignOp));
    
    /* 
     * insert "char** argv;" in globalScope,
     * and "globalArgv = argv;" in main function 
     */
    mainDef = mainDecl->get_definition();
    varDecl = buildVariableDeclaration(
                argvName,
                buildPointerType(buildPointerType(buildCharType())), 
                NULL,
                globalScope);
    insertStatement(getFirstStatement(globalScope), varDecl);
    assignOp = buildAssignOp(
            buildVarRefExp(argvName, globalScope),
            buildVarRefExp("argv", mainDef)); 
    mainDef->prepend_statement(buildExprStatement(assignOp));
}

/* use physis to represent stencils */
void autoModification(SgProject* project, StencilInfoList &stencilInfoList) {
    SgScopeStatement* globalScope = getFirstGlobalScope(project);
    dealWithArgcArgvInMain("globalArgc", "globalArgv", globalScope);

    dbg << GAP << "MODIFY stencils\n" << GAP; 

    /* add #include "physis/physis.h" */
	// myInsertHeader(project);

    /* build SgClassTypeMap */
    SgTypedefTypeMap typedefMap;
    buildSgTypedefTypeMap(typedefMap, project);
   
    int stencilId = 0;
    StrList idList;
    idList.push_back("i");
    idList.push_back("j");
    idList.push_back("k");
    FOREACH(StencilInfoList, stencilInfoList, iter) {
        modifyStencil(&(*iter), typedefMap, ++stencilId, idList, globalScope);
    } 
}

/* preserve the stencils of specified dimension only */
void checkDimension(StencilInfoList &stencilInfoList, int dimension) {
    if (dimension <= 0) return;

    FOREACH(StencilInfoList, stencilInfoList, iter) {
        if (iter->st_dim != dimension)
            stencilInfoList.erase(iter);
    }
}

} /* namespace physis */

static const char help_string[] = 
"Options:"
"\n     -m --modify     : Modify the stencil using DSL"
"\n     -d --dimension  : The dimension of the stencil to be modified, useful only when -s is set"
"\n     -h --help       : Print this help"
"\n     -I --include    : Pass include director to ROSE compiler"
"\n     -D --define     : Pass #define to ROSE compiler";

static const char options[] = "md:hI:D:";

static struct option opt_options[] = {
    {"modify", no_argument, 0, 'm'},
    {"dimension", required_argument, 0, 'd'},
    {"help", no_argument, 0, 'h'},
    {"include", required_argument, 0, 'I'},
    {"define", required_argument, 0, 'D'}
};

using namespace physis;
int main (int argc, char** argv)
{
    bool modify = false;
    int dimension = -1;

    /* solve the options */
    while (true) {
        int sscanf_return;
        int optchar = getopt_long(argc, argv, options, opt_options, NULL);
        if (optchar == -1)
            break;
        switch (optchar) {
            case 'm':
                modify = true;
                break;
            case 'd':
                sscanf_return = sscanf(optarg, "%d", &dimension);
                if (sscanf_return == EOF || sscanf_return == 0 
                        || dimension <= 0) {
                    fprintf(stderr, "Please enter a positive integer number for the dimesion of the stencil to be modified instead of \"-%c %s\"\n", optchar, optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 'h':
                printf("Usage: %s <options>\n%s\n", argv[0], help_string);
                return EXIT_SUCCESS;
            default:
                printf("%s", optarg);
                break;
        }
    }

    // Build the AST used by ROSE
    SgProject* project = frontend(argc, argv);
    project->skipfinalCompileStep(true);
    
    // Insert your own manipulations of the AST here...		

	/* get ALL for statements */
	SgNodeList forStmtList;
	getForStmtList(forStmtList, project); 
	
	/* get OUTEST for statements */
	SgNodeList outestForList;
	getOutestForList(outestForList, forStmtList);

    /* find STENCILS */
    SgNodeList stencilNodeList;
    getStencilNodeList(stencilNodeList, outestForList);
  
    /* get field length of struct */
    StructFieldMap structFieldMap;
    fillStructFieldLength(project, structFieldMap);

    /* get stencil INFO */
    StencilInfoList stencilInfoList;
    getStencilInfoList(stencilInfoList, stencilNodeList, structFieldMap);

    /* MODIFY stencils */
    if (modify) {
        checkDimension(stencilInfoList, dimension);
        autoModification(project, stencilInfoList);
    }

    // Run internal consistency tests on AST
    AstTests::runAllTests(project);
    
    // Generate source code from AST and invoke your
    // desired backend compiler
    return backend(project);
}

