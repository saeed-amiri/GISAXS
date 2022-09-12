import numpy as np
import pandas as pd


class Doc:
    """Calculate the reciprocal points of the structure
    Input:
        LAMMPS data file read by read_lmp_data
    Output:
        pd.DataFrame of the points reciprocal space
    """


class MakeReciprocal:
    """Mathematical caluatations of the space transformation"""
    def __init__(self, df: pd.DataFrame) -> None:
        reals_arr: np.array = self.df_to_array(df)  # Coords in real space
        self.to_reciprocal(reals_arr)
        del df

    def df_to_array(self, df: pd.DataFrame) -> np.array:
        """convert df to numpy array to calculate data file"""
        df_tmp: pd.DataFrame  # To save the coordinates and indices
        df_tmp = df[['atom_id', 'x', 'y', 'z']].copy()
        del df
        tmp_arr: np.array  # To convert, temporary
        tmp_arr = df_tmp.to_numpy(dtype=float)
        rows: int  # Number of rows in the Data file = NAtoms
        cols: int  # Number of columns at each row, then subtract id
        reals_arr: np.array  # To save the coords based on the id
        rows, cols = tmp_arr.shape
        reals_arr = np.zeros((rows, cols-1), dtype=float)
        for row in range(rows):
            id = int(tmp_arr[row][0]) - 1
            reals_arr[id][0] = tmp_arr[row][1]  # x positions
            reals_arr[id][1] = tmp_arr[row][2]  # y positions
            reals_arr[id][2] = tmp_arr[row][3]  # z positions
        del tmp_arr
        return reals_arr

    def to_reciprocal(self, reals_arr: np.array) -> None:
        """transfrom data to the reciprocal space"""
        pass
