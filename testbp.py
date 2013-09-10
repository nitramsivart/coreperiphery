import decellebeliefprop as bp
import numpy as np

def test_update_parameters():
  connected_adj = [[0,1,1],[1,0,1],[1,1,0]]
  gamma = [.5,.5]
  omega = [[1,0],[0,1]]
  message_field = ((1,0),
                   (0,1),
                   (0,1)) 
  messages = (
              (
                (1,0),
                (0,1),
                (0,1)
              ),
              (
                (1,0),
                (0,1),
                (0,1)
              ),
              (
                (1,0),
                (0,1),
                (0,1)
              ),
             )


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
  print bp.get_message_field(connected_adj, messages, omega, gamma)
  print bp.update_parameters(connected_adj, messages, message_field, joint_dist, omega)
  print bp.update_parameters2(connected_adj, messages, message_field, joint_dist, omega)

test_update_parameters()
