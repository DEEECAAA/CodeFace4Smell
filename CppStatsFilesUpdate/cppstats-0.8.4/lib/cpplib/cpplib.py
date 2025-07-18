'''
Created on Jun 11, 2010

@author: joliebig
'''

import sys

# pyparsing module
try:
    pversion = sys.version_info[0]

    if pversion == 2: import lib.pyparsing.pyparsing_py2 as pypa
    else: import lib.pyparsing.pyparsing_py3 as pypa
    pypa.ParserElement.enablePackrat()        # speed up parsing
    sys.setrecursionlimit(2000)               # handle larger expressions
except ImportError:
    print("pyparsing module not found! (python-pyparsing)")
    print("see http://pyparsing.wikispaces.com/")
    print("programm terminating ...!")
    sys.exit(-1)

# possible operands:
#   - hexadecimal number
#   - decimal number
#   - identifier
#   - macro function, which is basically expanded via #define
#     to an expression
__string = \
        pypa.Literal('\'').suppress() + \
        pypa.Word(pypa.alphanums+'_\\') + \
        pypa.Literal('\'').suppress()

__hexadec = \
        pypa.Literal('0x').suppress() + \
        pypa.Word(pypa.hexnums).\
        setParseAction(lambda t: str(int(t[0], 16))) + \
        pypa.Optional(pypa.Suppress(pypa.Literal('L')))

__integer = \
        pypa.Optional('~') + \
        pypa.Word(pypa.nums+'-').setParseAction(lambda t: str(int(t[0]))) + \
        pypa.Optional(pypa.Suppress(pypa.Literal('U'))) + \
        pypa.Optional(pypa.Suppress(pypa.Literal('L'))) + \
        pypa.Optional(pypa.Suppress(pypa.Literal('L')))

__identifier = \
        pypa.Word(pypa.alphanums+'_'+'-')
__arg = pypa.Word(pypa.alphanums+'_')
__args = __arg + pypa.ZeroOrMore(pypa.Literal(',').suppress() + \
        __arg)
__fname = pypa.Word(pypa.alphas, pypa.alphanums + '_')
__function = pypa.Group(__fname + pypa.Literal('(').suppress() + \
        __args + pypa.Literal(')').suppress())


def _parseIfDefExpression(ifdefexp):
    """This function parses a given ifdef-expression and
    rewrites the expression according to the given __pt mapping.
    This one is used to make use of a csp solver without using
    a predicate."""
    mal = list()

    def _rewriteOne(param):
        """This function returns each one parameter function
        representation for csp."""
        op, ma = param[0]
        mal.append(ma)
        if op == '!': ret = op + '(' + ma + ')'
        if op == 'defined': ret = ma
        return  ret

    def _rewriteTwo(param):
        """This function returns each two parameter function
        representation for csp."""
        mal.extend(param[0][0::2])
        ret = param[0][1]
        ret = '(' + ret.join(map(str, param[0][0::2])) + ')'
        return ret

    operand = __string | __hexadec | __integer | \
            __function | __identifier
    operators = pypa.oneOf('&& ||') # extend with furhter operators
    expr = pypa.operatorPrecedence(operand, [
        ('defined', 1, pypa.opAssoc.RIGHT, _rewriteOne),
        ('!',  1, pypa.opAssoc.RIGHT, _rewriteOne),
        (operators, 2, pypa.opAssoc.LEFT, _rewriteTwo),
    ]) + pypa.StringEnd()

    try:
        rsig = expr.parseString(ifdefexp)[0]
    except pypa.ParseException as e:
        print(('ERROR (parse): cannot parse sig (%s) -- (%s)' %
                (ifdefexp, e.col)))
        return ifdefexp
    except RuntimeError:
        print(('ERROR (time): cannot parse sig (%s)' % (ifdefexp)))
        return ifdefexp
    return (mal, ''.join(rsig))


__ifdefexplist = []
def _collectIfdefExpressions(fname):
    '''
    This method filters all ifdef expressions out of a file and returns them as a list.
    '''

    def _extractIfdefExpression(tokens):
        global __ifdefexplist
        __ifdefexplist += tokens

    __macro = pypa.Literal('#') \
            + pypa.oneOf("if ifdef ifndef elif") \
            + pypa.OneOrMore(pypa.Word(pypa.alphanums+"&|><=^")) \
                    .setParseAction(_extractIfdefExpression) \
            + pypa.StringEnd()

    with open(fname, 'r') as fd:
        for line in fd:
            try:
                print((__macro.parseString(line)))
            except pypa.ParseException:
                pass
    return __ifdefexplist

def _filterAnnotatedIfdefs(fnamein, fnameout):
    '''
    This method removes all preprocessor annotated lines from the input.
    '''
    inifdef = 0

    with open(fnameout, 'w') as fdout:
        with open(fnamein, 'r') as fdin:
            for line in fdin.readlines():
                # line is a preprocessor directive; determine weather its a conditional inclusion
                # directive or something else
                if line.startswith('#'):
                    parseddirective = line.split(' ', 1)
                    directive = parseddirective[0].strip()
                    if directive in ['#if', '#ifdef', '#ifndef']:
                        inifdef += 1
                    elif directive in ['#else', '#elif']:
                        pass
                    elif directive == '#endif':
                        inifdef -= 1
                    elif directive in ['#line', '#error', '#pragma', '#define', '#undef', '#include', '#ident', '#warning', '#include_next']:
                        if inifdef:
                            fdout.write('\n')
                            continue
                    else:
                        print(("ERROR: directive (%s) not processed!" % parseddirective))
                    fdout.write(line)
                # found regular C code
                else:
                    fdout.write(line)


##################################################
if __name__ == '__main__':
    symbols, expression = _parseIfDefExpression('AA && BB')
    print((symbols, expression))
    print((_collectIfdefExpressions('/home/joliebig/workspace/reverse_cpp/test/test.c')))
    _filterAnnotatedIfdefs('/home/joliebig/workspace/reverse_cpp/test/filterannotateddirectives.h',
                           '/home/joliebig/workspace/reverse_cpp/test/filterannotateddirectives_out.h')
