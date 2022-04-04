\ Model matlab
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  x
Subject To
 R0: x + y + z = 1
 std_cone: [ x ^2 + y ^2 - z ^2 ] <= 0
 rot_cone: [ x ^2 - y * z ] <= 0
Bounds
End
