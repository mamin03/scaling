from __future__ import division
from scitbx.lstbx import normal_eqns
from scitbx.lstbx import normal_eqns_solving
import boost.python
ext = boost.python.import_ext("xscale_ext")
from xscale_ext import *

