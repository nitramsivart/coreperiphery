import oldbeliefprop as bp
import numpy as np

def test_update_parameters():
  connected_adj = [[0,1,1],[1,0,1],[1,1,0]]
  gamma = [.5,.5]
  omega = [[1,0],[0,1]]
  message_field = ((1,0),
                   (0,1),
                   (0,1)) 
  messages = np.zeros([3,3,2])
  joint_dist = (
                 (
                   None,
                   (
                    (0,1),
                    (0,0)
                   ),
                   (
                    (0,1),
                    (0,0)
                   )
                 ),
                 (
                   (
                     (0,0),
                     (1,0)
                   ),
                   None,
                   (
                     (0,0),
                     (0,1)
                   )
                 ),
                 (
                   (
                     (0,0),
                     (1,0)
                   ),
                   (
                     (0,0),
                     (0,1)
                   ),
                   None
                 )
               )
  print bp.update_parameters(connected_adj, messages, message_field, joint_dist)

test_update_parameters()
