# allow a class based construction for weave files
import numpy as np

def get_ctype(vname,loc,glob):
  base_type = eval('type({})'.format(vname),glob,loc)
  if base_type == np.int:
    c_type = 'int'
  elif base_type in [np.float64,np.float]:
    c_type = 'double'
  elif base_type == np.ndarray:
    np_type = eval('{}.dtype'.format(vname),glob,loc)
    if np_type == np.int:
      c_type = 'long int*'
    elif np_type == np.float:
      c_type = 'double*'
    else:
      print 'Unknown numpy type {}'.format(np_type)
      c_type = None
  else:
    print 'Unknown python type {}'.format(base_type)
    c_type = None
  return c_type

def generate_code(fname,args,loc,glob,class_name='WeaveContainerClass',dyn_flag='XXX_DYNAMIC_CODE_XXX',run='run',inst_name='wcc_inst',store=None):
  var_types = [get_ctype(iarg,loc,glob) for iarg in args]
  var_assign = ''.join(['this->{} = {};\n'.format(iarg,iarg) for iarg in args])
  var_declare = [(itype+' '+iarg) for (itype,iarg) in zip(var_types,args)]
  var_constructor = class_name + '(' + ', '.join(var_declare) + ') {\n' + var_assign + '}\n\n'
  var_define = ''.join([iarg+';\n' for iarg in var_declare])

  dyn_setup = var_constructor + var_define
  dyn_execute = class_name + ' ' + inst_name + ' = ' + class_name + '(' + ','.join(args) + ');\n' + inst_name + '.' + run + '();\n'
  base_code = open(fname,'r').read()
  code = base_code.replace(dyn_flag,dyn_setup) + '\n' + dyn_execute

  if store:
    code_store = open(store,'w+')
    code_store.write(code)
    code_store.close()

  return code
