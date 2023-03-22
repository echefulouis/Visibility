import matplotlib.pyplot as plt

class ObstaclePolygon:

    def __init__(self, x_list, y_list):
        self.x_list = x_list
        self.y_list = y_list

        self.close_polygon()
        self.make_clockwise()

    def make_clockwise(self):
        if not self.is_clockwise():
            self.x_list = list(reversed(self.x_list))
            self.y_list = list(reversed(self.y_list))

    def is_clockwise(self):
        n_data = len(self.x_list)
        eval_sum = sum([(self.x_list[i + 1] - self.x_list[i]) *
                        (self.y_list[i + 1] + self.y_list[i])
                        for i in range(n_data - 1)])
        eval_sum += (self.x_list[0] - self.x_list[n_data - 1]) * \
                    (self.y_list[0] + self.y_list[n_data - 1])
        return eval_sum >= 0

    def close_polygon(self):
        is_x_same = self.x_list[0] == self.x_list[-1]
        is_y_same = self.y_list[0] == self.y_list[-1]
        if is_x_same and is_y_same:
            return  # no need to close

        self.x_list.append(self.x_list[0])
        self.y_list.append(self.y_list[0])

    def plot(self):
        plt.plot(self.x_list, self.y_list, "-k")