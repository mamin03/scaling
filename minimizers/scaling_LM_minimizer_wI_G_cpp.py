from __future__ import division
from scitbx.array_family import flex
from cxi_xdr_xes.XScale import xscale_G_wI_LM
from scitbx.lstbx import normal_eqns

class xscale(xscale_G_wI_LM, normal_eqns.non_linear_ls_mixin):

  def __init__(self, I, w, hkl, frames, G, Ih):
    self.I = flex.double(I)
    self.w = flex.double(w)
    self.hkl = flex.size_t(hkl)
    self.frames = flex.size_t(frames)
    self.G = G
    self.Ih = flex.double(Ih)
    self.S = flex.double(Ih)
    self.n_obs= len(I)
    self.n_frames = len(G)
    self.n_prm = self.n_frames + len(self.S)
    print "Number of parameters :", self.n_prm
    self.count=0
    super(xscale, self).__init__(self.n_prm)
    self.set_cpp_data(len(self.I))
    self.restart()
    self.build_up()

  def restart(self):
    self.x = flex.double(self.n_prm,1)
    self.old_x = None
    print "Restart"

  def parameter_vector_norm(self):
    return self.x.norm()

  def build_up(self, objective_only=False):
    print "Trial number ...", self.count
    self.reset()
    Gm = self.x[:self.n_frames]
    Ih = self.x[self.n_frames:]
    self.build_up_normal_matrix(objective_only, self.I,
                                self.w,
                                self.hkl,
                                self.frames,
                                self.S,
                                Gm,
                                Ih,
                                self.n_frames)
#   self.build_up_sparse_normal_matrix(objective_only, self.I,
#                                      self.w,
#                                      self.hkl,
#                                      self.frames,
#                                      self.S,
#                                      Gm,
#                                      Ih,
#                                      self.n_frames)
    self.count+=1
  def step_forward(self):
    self.old_x = self.x.deep_copy()
    self.x += self.step()

  def step_backward(self):
    assert self.old_x is not None
    self.x, self.old_x = self.old_x, None
