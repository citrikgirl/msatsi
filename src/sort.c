//-------------------------------------------------------------------------------------------------
void sort(float x[], int n) {
  int i;
  double hold;
  int flag;

  flag = 1;
  while (flag) {
    flag = 0;
    for (i = 0; i < n - 1; ++i) {
      if (x[i] <= x[i + 1])
        continue;
      flag = 1;
      hold = x[i];
      x[i] = x[i + 1];
      x[i + 1] = hold;
    }
  }
  return;
}
