class Read(object):
      def __init__(self, slako_name):

          self.slako_name = slako_name

      def read(self):

          fileIn = file(self.slako_name,"r")
          self.lines = fileIn.readlines()

          fileIn.close()

      def get_numbers(self):
 
          count = 0
          for line in self.lines:
              pieces = line.split()
              size = len(pieces)
              count = 0
              for ipieces in pieces:
                  count += 1
                  #
                  if ipieces == "20*0.0":
#                     print [0.0 for i in xrange(20)]
                     print ' '.join(["0.0" for i in xrange(20)])
                     break
                  #
                  elif ipieces == "19*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(19)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(19)])
                  #
                  elif ipieces == "18*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(18)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(18)])
                  #
                  elif ipieces == "17*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(17)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(17)])
                  #
                  elif ipieces == "16*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(16)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(16)])
                  #
                  elif ipieces == "15*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(15)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(15)])
                  #
                  elif ipieces == "14*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(14)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(14)])
                  #
                  elif ipieces == "13*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(13)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(13)])
                  #
                  elif ipieces == "12*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(12)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(12)])
                  #
                  elif ipieces == "11*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(11)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(11)])
                  #
                  elif ipieces == "10*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(10)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(10)])
                  #
                  elif ipieces == "9*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(9)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(9)])
                  #
                  elif ipieces == "8*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(8)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(8)])
                  # 
                  elif ipieces == "7*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(7)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(7)])
                  #
                  elif ipieces == "6*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(6)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(6)])
                  #
                  elif ipieces == "5*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(5)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(5)])
                  #
                  elif ipieces == "4*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(4)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(4)])
                  #
                  elif ipieces == "3*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(3)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(3)])
                  #
                  elif ipieces == "2*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(2)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(2)])
                  #
                  elif ipieces == "1*0.0":
                       if count != size:
                          print ' '.join(["0.0" for i in xrange(1)]),
                       if count == size:
                          print ' '.join(["0.0" for i in xrange(1)])
                  #
                  else:
                    if count != size:
                       print ipieces,
                    if count == size:
                       print ipieces
#----------------------------------------------
x = Read("MgMg")
x.read()
x.get_numbers()
