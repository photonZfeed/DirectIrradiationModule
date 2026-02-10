import numpy as np

class Grid:
    """
    Represents a two-dimensional, symmetric grid of points centered at the origin (0, 0).

    This class allows for the creation of a customizable 2D grid, where the user can specify the total width and height, the margin (space) to leave on the sides and top/bottom, and the step size (distance) between adjacent grid points. The grid points are generated as a NumPy array of (x, y) coordinates, suitable for geometric modeling, simulation, or spatial analysis tasks.

    **Usage Notes:**
        - The grid is always centered at (0, 0).
        - Margins are excluded from the grid area; no points are placed within the specified margins.
        - The step size determines the spacing between adjacent points and must be positive.
        - All dimensions and spacings are in arbitrary but consistent units (e.g., meters).
    """

    def __init__(self, width: float = 10.0, height: float = 10.0, side_space: float = 1.0, top_bottom_space: float = 1.0, step: float = 2.0):
        """
        Initialize a symmetric 2D grid with specified dimensions, margins, and spacing.

        :param width: Total width of the grid (units: arbitrary, e.g., meters). Must be positive. Default is 10.0.
        :type width: float
        :param height: Total height of the grid (units: arbitrary, e.g., meters). Must be positive. Default is 10.0.
        :type height: float
        :param side_space: Margin to leave on both left and right sides (units: same as width). Must be non-negative. Default is 1.0.
        :type side_space: float
        :param top_bottom_space: Margin to leave on both top and bottom sides (units: same as height). Must be non-negative. Default is 1.0.
        :type top_bottom_space: float
        :param step: Distance between adjacent grid points (units: same as width/height). Must be positive. Default is 2.0.
        :type step: float

        :raises ValueError: If any of the dimensions or step are not positive, or if margins are negative or too large for the grid size.

        **Side Effects:**
            - Initializes the `grid_points` attribute as a NumPy array of shape (N, 2).

        **Notes:**
            - Margins are subtracted from both sides (i.e., total margin is 2 * side_space for width, 2 * top_bottom_space for height).
            - The grid is always centered at (0, 0).
        """

        if width <= 0 or height <= 0:
            raise ValueError("Grid width and height must be positive.")
        if step <= 0:
            raise ValueError("Grid step size must be positive.")
        if side_space < 0 or top_bottom_space < 0:
            raise ValueError("Grid margins must be non-negative.")
        if 2 * side_space >= width:
            raise ValueError("Side margins are too large for the specified width.")
        if 2 * top_bottom_space >= height:
            raise ValueError("Top/bottom margins are too large for the specified height.")

        self.width = width
        self.height = height
        self.side_margin = side_space
        self.top_bottom_margin = top_bottom_space
        self.step = step
        self.grid_points = self.create_grid()

    def create_grid(self):
        """
        Generate the coordinates of all grid points within the specified dimensions and margins.

        :return: NumPy array of shape (N, 2), where each row is an (x, y) coordinate of a grid point.
        :rtype: numpy.ndarray

        **Notes:**
            - The grid is symmetric and centered at (0, 0).
            - Points are spaced by `step` units in both x and y directions.
            - Margins are respected; no points are placed within the specified margins.

        **Assumptions:**
            - All class attributes are valid and checked during initialization.
        """
        n_x = int((self.width / 2 - self.side_margin) // self.step)  # number of grid points in x direction (one side)
        n_y = int((self.height / 2 - self.top_bottom_margin) // self.step)  # number of grid points in y direction (one side)

        # Create grid points
        x = np.arange(-n_x * self.step, (n_x + 1) * self.step, self.step)
        y = np.arange(-n_y * self.step, (n_y + 1) * self.step, self.step)
        X, Y = np.meshgrid(x, y)
        grid_points = np.stack([X.ravel(), Y.ravel()], axis=1)

        return grid_points