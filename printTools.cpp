
void  printMat(vector<vector<double>>& mat) {
    int row = mat[0].size();
    int col = mat[1].size();
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            cout << mat[i][j] << " ";

        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------" << endl;
}

void  printMat(vector<vector<double>>& mat, string nameVar) {
    cout << "------------------------" << nameVar << "------------------------" << endl;

    int row = mat[0].size();
    int col = mat[1].size();
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            cout << mat[i][j] << " ";

        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------" << endl;
}

void  printVec(vector<double>& vec) {

    int row = vec.size();
    for (int i = 0; i < row; i++) {
        cout << vec[i] << endl;
    }
    cout << "---------------------------------------------------------------" << endl;
}
