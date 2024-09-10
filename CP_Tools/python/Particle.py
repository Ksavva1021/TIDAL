class Particle:
   """
   Description: Represents a particle with a Fourvector and electric charge.
   Attributes:
   - Fourvector (tuple): A four-dimensional vector describing the particle's properties.
   - charge (float): The electric charge of the particle.
   """
   def __init__(self, Fourvector, charge):
      self.Fourvector = Fourvector
      self.charge = charge