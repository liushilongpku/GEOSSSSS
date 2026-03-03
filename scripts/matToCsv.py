import numpy as np
import sys

def parse_input(input_file):
    data = []
    with open(input_file, 'r') as f:
        # 跳过表头
        next(f)
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    global_row_id = int(parts[1])
                    global_col_id = int(parts[2])
                    value = float(parts[3])
                    data.append((global_row_id, global_col_id, value))
                except ValueError:
                    continue
    return data

def write_to_csv(matrix, output_file):
    """使用numpy写入CSV文件"""
    np.savetxt(output_file, matrix, delimiter=',', fmt='%.6f')

def main(input_file, output_file):
    # 解析输入数据
    data = parse_input(input_file)
    
    if not data:
        print("No valid data found in input file.")
        return
    
    # 计算最大行ID和列ID
    max_row = max(row for row, col, val in data)
    max_col = max(col for row, col, val in data)
    
    # 创建零矩阵
    matrix = np.zeros((max_row + 1, max_col + 1))
    
    # 填充数据
    for row, col, val in data:
        matrix[row, col] = val
    
    # 根据输出文件扩展名选择格式
    if output_file.endswith('.xlsx'):
        try:
            import pandas as pd
            df = pd.DataFrame(matrix)
            df.to_excel(output_file, index=False, header=False)
            print(f"Data successfully written to {output_file}")
        except ImportError:
            # 如果没有pandas，保存为CSV
            csv_file = output_file.replace('.xlsx', '.csv')
            write_to_csv(matrix, csv_file)
            print(f"pandas not found. Data written to {csv_file} instead")
    elif output_file.endswith('.csv'):
        write_to_csv(matrix, output_file)
        print(f"Data successfully written to {output_file}")
    else:
        # 默认保存为CSV
        write_to_csv(matrix, output_file)
        print(f"Data successfully written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python matToCsv.py <input_file> <output_file>")
    else:
        main(sys.argv[1], sys.argv[2])
