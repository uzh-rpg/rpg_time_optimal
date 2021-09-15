import yaml

class Track:
  def __init__(self, filename = ""):
    self.init_pos = None
    self.init_att = None
    self.init_vel = None
    self.init_omega = None
    self.end_pos = None
    self.end_att = None
    self.end_vel = None
    self.end_omega = None
    self.gates = []
    self.ring = False
    if filename:
      self.load(filename)

  def addGate(self, gate):
    self.gates += [gate]
    print(self.gates)

  def load(self, filename):
    print("Loading track from " + filename)
    with open(filename, 'r') as file:
      track = yaml.load(file, Loader=yaml.FullLoader)
      
    if 'gates' in track:
      self.gates = track['gates']
    else:
      print("No gates specified in " + filename)
    if 'initial' in track:
      initial = track['initial']
    else:
      initial = track
    if 'position' in initial:
      self.init_pos = initial['position']
    if 'attitude' in initial:
      self.init_att = initial['attitude']
    if 'velocity' in initial:
      self.init_vel = initial['velocity']
    if 'omega' in initial:
      self.init_omega = initial['omega']
    if 'end' in track:
      end = track['end']
      if 'position' in end:
        self.end_pos = end['position']
      if 'attitude' in end:
        self.end_att = end['attitude']
      if 'velocity' in end:
        self.end_vel = end['velocity']
      if 'omega' in end:
        self.end_omega = end['omega']
    if 'ring' in track:
      self.ring = track['ring']