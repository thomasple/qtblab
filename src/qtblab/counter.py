class Counter:
  def __init__(self,nseg,startsave=1):
    self.__reset_step=True
    self.__i=0
    self.__i_avg=0
    self.__startsave=startsave
    self.__nseg=nseg

  def i(self):
    return self.__i
  
  def nsample(self):
    return max(self.__i_avg-self.__startsave+1,1)

  def is_reset_step(self):
    return self.__i == 0

  def increment(self):
    self.__i+=1
    if(self.__i>=self.__nseg):
      self.__i=0
      self.__i_avg+=1