import sys
from PyQt5 import QtGui, QtCore, QtWidgets


class MyWidget(QtWidgets.QWidget):
    def __init__(self, title=None, savefilename='mydrawingx.jpg', spiralx=[], labels=[], n=0):
        super(MyWidget, self).__init__()
        self.title = title
        self.spiral = spiralx
        self.labels = labels
        self.n = n
        if '>' in savefilename:
            self.fname = savefilename[1:]
        else:
            self.fname = savefilename

        self.initUI()
        self.test = self
        self.label = QtWidgets.QLabel(self)

    def initUI(self):
        #self.setGeometry(100, 400, 600, 400)
        self.setGeometry(300, 300, 650, 450)
        self.setWindowTitle(self.title)
        self.show()

    def paintEvent(self, e):
        qp = QtGui.QPainter()
        self.drawspiral(qp)




    def drawRectes(self, qp):
        '''
        # tonos de azul
        qp.setBrush(QtGui.QColor(0, 255, 255))
        qp.drawRect(0, 310, 40, 40)
        # Blue
        qp.setBrush(QtGui.QColor(3, 169, 244)) # blue
        qp.drawRect(40, 310, 80, 40)
        # red
        qp.setBrush(QtGui.QColor(3,169, 244))
        qp.drawRect(80, 310, 40, 40)
        '''

        # Red
        qp.setBrush(QtGui.QColor(221, 34, 0))
        qp.drawRect(10, 15, 90, 60)

        # Yellow
        qp.setBrush(QtGui.QColor(255, 193, 7))
        qp.drawRect(110, 15, 90, 60)

        # Gray
        qp.setBrush(QtGui.QColor(158, 158, 158))
        qp.drawRect(10, 85, 90, 60)

        # Green
        qp.setBrush(QtGui.QColor(76, 175, 80))
        qp.drawRect(110, 85, 90, 60)

        # LightGreen
        qp.setBrush(QtGui.QColor(0, 255, 0))
        qp.drawRect(10, 155, 90, 60)

        qp.setBrush(QtGui.QColor(3, 169, 244))  # blue
        qp.drawRect(110, 155, 90, 60)

    def drawspiral(self, qp):
        offsetx = 5
        offsety = 5
        sizex = 50
        sizey = 50
        multiplicator = len(self.spiral)

        w, h = multiplicator * sizex + 2 * offsetx, multiplicator * sizey + 2 * offsety
        self.resize(w, h)
        pixmap = QtGui.QPixmap(w, h)
        pixmap.fill(QtCore.Qt.white)
        qp.begin( pixmap)
        qp.setBrush(QtGui.QColor(255, 255, 255))  # white
        qp.drawRect(300, 300, w, h)
         # fill colour
        for i in range(multiplicator):
            for j in range(multiplicator):
                if self.spiral[i][j] == 0:
                    qp.setBrush(QtGui.QColor(255, 255, 255))
                    qp.setPen(QtGui.QColor(255, 255, 255))
                elif self.spiral[i][j] == 1:
                    qp.setBrush(QtGui.QColor(3, 169, 244))
                    qp.setPen(QtGui.QColor(3, 169, 244))
                elif self.spiral[i][j] == 2:
                    qp.setBrush(QtGui.QColor(221, 34, 0))
                    qp.setPen(QtGui.QColor(221, 34, 0))
                elif self.spiral[i][j] == 3:
                     qp.setBrush(QtGui.QColor(255, 193, 7))
                     qp.setPen(QtGui.QColor(255, 193, 7))
                elif self.spiral[i][j] == 4:
                    qp.setBrush(QtGui.QColor(0, 255, 0))
                    qp.setPen(QtGui.QColor(0, 255, 0))
                elif self.spiral[i][j] == 5:
                    qp.setBrush(QtGui.QColor(158, 158, 158))
                    qp.setPen(QtGui.QColor(158, 158, 158))

                qp.drawRect(j*sizex+offsetx, i*sizey+offsety, sizex,sizey)

                if self.labels[i][j] == '--':
                    self.labels[i][j] = ''

                try:
                    qp.setPen(QtGui.QColor(255, 255, 255))
                    p = QtCore.QPointF()
                    p.setX(j * sizex + 2 * offsetx + 5)
                    p.setY(i * sizey + 3 * offsety + 20)
                    qp.drawText(p, self.labels[i][j])
                except:
                    pass

        # Spiral line over rectangle edge
        qp.setBrush(QtGui.QColor(255, 255, 255))
        qp.setPen(QtGui.QColor(255, 255, 255))
        pen = QtGui.QPen()
        pen.setColor(QtGui.QColor(255, 255, 255))
        pen.setWidth(5)
        qp.setPen(pen)
        four_states_machine = ["Down", "left", "up", "right"]
        movement = [(0, 1), (-1, 0), (0, -1), (1, 0)]
        stepf = 0
        longitud = 1
        posrefx, posrefy = w / 2, h / 2
        if multiplicator % 2 != 0:
            posrefx, posrefy = posrefx + sizex/2, posrefy-sizey/2

        while stepf < 2 * multiplicator:
            x = movement[stepf % 4]
            newpx = posrefx + sizex * x[0] * longitud
            newpy = posrefy + sizey * x[1] * longitud
            qp.drawLine(posrefx, posrefy, newpx, newpy)
            posrefx = newpx
            posrefy = newpy
            stepf += 1
            x = movement[stepf % 4]
            newpx = posrefx + sizex * x[0] * longitud
            newpy = posrefy + sizey * x[1] * longitud
            qp.drawLine(posrefx, posrefy, newpx, newpy)
            posrefx = newpx
            posrefy = newpy
            stepf += 1
            longitud += 1
        qp.end()
        self.label.setPixmap(pixmap)
        self.label.show()
        self.label.pixmap().save(self.fname)


def main():
    myspiral = [[1, 0, 1, 0, 1, 0],
                [0, 2, 0, 2, 0, 2],
                [3, 0, 3, 0, 3, 0],
                [0, 4, 0, 4, 0, 4],
                [1, 0, 2, 0, 3, 0],
                [0, 1, 0, 4, 0, 1]]

    mylabels = [['T26', '27-', 'T28', '29-', 'T30', '31-'],
                ['25-', 'G10', '11-', 'G12', '13-', 'G32'],
                ['A24', '9-', 'A2', '3-', 'A14', '33-'],
                ['23-', 'C8', '1-', 'C4', '15-', 'C34'],
                ['T22', '7-', 'G6', '5-', 'A16', '35-'],
                ['21-', 'T20', '19-', 'C18', '17-', 'T36'], ]

    app = QtWidgets.QApplication(sys.argv)
    ex = MyWidget(title='Four colors map representation', spiralx=myspiral, labels=mylabels)

    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
