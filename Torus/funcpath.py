import operator
from decimal import Decimal
import numpy as np
import plotly.graph_objects as go

class FuncPath:
    """
    A path between two points on the Torus, to which the functional iteration is applied.


    Attributes
    ----------
    u : numpy.ndarray
        Array of values between 0 and 2pi in toroidal direction

    v: numpy.ndarray
        Array of values between 0 and 2pi in poloidal direction

    period: float
        Value 2pi. Period of Torus parameters

    length: int
        length of u

    current_iteration:int
        how many iteractions has this path gone through since the initial guess

    c : float
    a : float
        coeffients of torus parametrization. e.g x = (c + a*cos(v)) * cos(u)

    


    Methods
    -------
    periodic_op(arr1, arr2, op_type)
        Element-wise operation on arrays with consideration of the period
    periodic_op_scalar(a,b,op_type)
        Scalar operation with consideration of the period
    periodic_caliberate(arr)
        Recalibrate array with respect to the period
    functional_iteration()
        Executes one round of functional iteration.
    curve_length()
        Calculate arc length of path
    tor2cart()
        Calculate cartesian coordinates from toroidal coordinates
    draw_path(ax, prev_path=None)
        Draw path on existing matplotlib figure     
    plotly_draw_path_go(fig)
        Draw path on existing Dash figure, graph_objects style
    plotly_draw_path(fig)
        Draw path on existing Dash figure, dictionary style
    plotly_update_path(i, fig)
        Update path on Dash figure, dictionary style
    """

    period = 2 * np.pi

    def __init__(self,u,v,c,a):
        self.u = u
        self.v = v
        self.c = c
        self.a = a
        self.length = len(u)
        self.N = self.length - 1 # number of intervals in the path
        self.current_iteration = 0

    def periodic_op(self,arr1, arr2, op_type):
        """Element-wise operation on arrays with consideration of the period.
        
        Parameters
        ----------
        arr1: numpy.ndarray
            Subslice of an array from the third element to the last element

        arr2: numpy.ndarray
            Subslice of the same array from the first element to the third last element.

        op_type: builtin_function_or_method
            A math operator function such as operator.add

        
        Returns
        -------
        result: numpy.ndarray
        """

        # TODO: does this yield all values with magnitude less than 2pi? does it need to ?

        result = np.copy(arr1)

        for i in range(arr1.size):
            a = arr1[i]
            b = arr2[i]
            if abs(b-a) > self.period / 2:
                if a < b:
                    a = a + self.period
                else:
                    b = b + self.period
            result[i] = op_type(a,b)

        return result


    def periodic_op_scalar(self,a,b,op_type):
        """ Scalar operation with consideration of the period
        
        Parameters
        ----------
        a: float
        b: float

        Returns
        -------
        float
        """
        if abs(b-a) > self.period / 2:
            if a < b:
                a = a + self.period
            else:
                b = b + self.period
        return op_type(a,b)


    def periodic_caliberate(self,arr):
        """Recalibrate array with respect to the period
        

        Parameters
        ----------
        arr: numpy.ndarray


        Returns
        -------
        result: numpy.ndarray
        """

        result = np.copy(arr)
        for i in range(arr.size):
            result[i] = float(Decimal(arr[i]) % Decimal(self.period))
        return result


    def functional_iteration(self):
        """Executes one round of functional iteration."""

        u_diff = self.periodic_op(self.u[2:], self.u[0:-2], operator.sub)
        v_diff = self.periodic_op(self.v[2:], self.v[0:-2], operator.sub)
        u_sum = self.periodic_op(self.u[2:], self.u[0:-2], operator.add)
        v_sum = self.periodic_op(self.v[2:], self.v[0:-2], operator.add)

        new_u = np.copy(self.u)
        frac1 = np.divide(self.a * np.sin(self.v), self.c + self.a * np.cos(self.v))
        new_u[1:-1] = u_sum/2 + np.multiply( frac1[1:-1], np.multiply( u_diff, v_diff )) / 4

        new_v = np.copy(self.v)
        frac2 = np.multiply(np.sin(self.v)/self.a, self.c + self.a*np.cos(self.v))
        new_v[1:-1] = v_sum/2 + np.multiply( frac2[1:-1] , np.multiply(u_diff,u_diff) ) / 8

        self.u = self.periodic_caliberate(new_u)
        self.v = self.periodic_caliberate(new_v)

        self.current_iteration += 1


    def curve_length(self):
        # the path is parametrized by arclength??
        """Calculate curve length of path
        
        Returns
        -------
        float
        """

        delta = 1/self.N
        delta2 = 2/self.N
        
        v_diff_0 = self.periodic_op_scalar(self.v[1], self.v[0], operator.sub)
        u_diff_0 = self.periodic_op_scalar(self.u[1], self.u[0], operator.sub)
        u_diff_end = self.periodic_op_scalar(self.u[-1], self.u[-2], operator.sub)
        v_diff_end = self.periodic_op_scalar(self.v[-1], self.v[-2], operator.sub)

        temp1 = np.sqrt( ( self.a * v_diff_0 / delta )**2 + ( (self.c + self.a*np.cos(self.v[0])) * u_diff_0 / delta )**2)
        temp2 = np.sqrt( ( self.a * v_diff_end / delta )**2 + ( (self.c + self.a*np.cos(self.v[0])) * u_diff_end / delta )**2)
        sum = (temp1 + temp2)/2

        for i in range(1,self.N):
            v_diff_i = self.periodic_op_scalar(self.v[i+1], self.v[i-1], operator.sub)
            u_diff_i = self.periodic_op_scalar(self.u[i+1], self.u[i-1], operator.sub)

            sum += np.sqrt( ( self.a * v_diff_i / delta2 )**2 + ( (self.c + self.a*np.cos(self.v[i])) * u_diff_i / delta2 )**2)

        return sum * delta


    def tor2cart(self):
        """Calculate cartesian coordinates from toroidal coordinates

        Returns
        -------
        x, y, z : Array of float
        """
        x = (self.c + self.a*np.cos(self.v)) * np.cos(self.u)
        y = (self.c + self.a*np.cos(self.v)) * np.sin(self.u)
        z = self.a * np.sin(self.v)
        return x, y, z

    def draw_path(self, ax, prev_path=None): # 
        """
        Draw path on existing matplotlib figure       
        """
        if prev_path is not None:
            prev_path.remove()

        x,y,z = self.tor2cart()
        return ax.plot(x, y, z)

    def plotly_draw_path_go(self,fig):
        """Draw path on existing Dash figure, graph_objects style"""
        x, y, z = self.tor2cart()
        fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', name='lines'))

    def plotly_draw_path(self, fig):
        """Draw path on existing Dash figure, dictionary style"""
        x, y, z = self.tor2cart()
        new = {'mode': 'lines', 'name': 'lines', 'x': x, 'y': y, 'z': z, 'type': 'scatter3d'}
        fig['data'].append(new)

    def plotly_update_path(self, i, fig):
        """Update path on Dash figure, dictionary style"""
        x, y, z = self.tor2cart()
        new = {'mode': 'lines', 'name': 'lines', 'x': x, 'y': y, 'z': z, 'type': 'scatter3d'}
        fig['data'][i+1] = new


