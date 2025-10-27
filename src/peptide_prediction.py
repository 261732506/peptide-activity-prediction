"""
多功能肽智能改造系统 - 优化版
基于ESM-2-3B模型的功能位点识别与序列编辑
优化内存使用，支持大规模批处理
"""

import torch
import esm
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.metrics import classification_report, f1_score
import warnings
warnings.filterwarnings('ignore')
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import gc
import os

# 设置环境变量优化内存分配
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'

class ESMPeptideAnalyzer:
    """基于ESM-2的肽功能分析器 - 优化内存版"""
    
    def __init__(self, device='cuda' if torch.cuda.is_available() else 'cpu', batch_size=32):
        """
        初始化ESM-2-3B模型
        batch_size: 根据96G显存优化的批大小
        """
        print(f"Loading ESM-2-3B model on {device}...")
        self.model, self.alphabet = esm.pretrained.esm2_t36_3B_UR50D()
        self.model = self.model.to(device)
        self.model.eval()
        self.batch_converter = self.alphabet.get_batch_converter()
        self.device = device
        self.batch_size = batch_size
        
        # 打印显存信息
        if device == 'cuda':
            print(f"GPU Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.2f} GB")
            print(f"Allocated: {torch.cuda.memory_allocated() / 1e9:.2f} GB")
            
    def get_embeddings_and_attention_batch(self, sequences: List[str]) -> Tuple[np.ndarray, List[np.ndarray]]:
        """
        批处理获取序列的嵌入和注意力权重
        返回：(embeddings, attention_weights)
        """
        all_embeddings = []
        all_attention_weights = []
        
        # 分批处理
        for i in tqdm(range(0, len(sequences), self.batch_size), desc="Processing batches"):
            batch_seqs = sequences[i:i + self.batch_size]
            
            # 准备批次数据
            data = [(f"seq_{j}", seq) for j, seq in enumerate(batch_seqs)]
            batch_labels, batch_strs, batch_tokens = self.batch_converter(data)
            
            # 限制最大序列长度以节省内存
            max_len = min(50, batch_tokens.shape[1])
            batch_tokens = batch_tokens[:, :max_len].to(self.device)
            
            # 提取特征
            with torch.no_grad():
                # 不返回contacts以节省内存
                results = self.model(batch_tokens, repr_layers=[36], return_contacts=False)
                
                # 提取嵌入（最后一层）
                token_embeddings = results["representations"][36]
                
                # 处理每个序列
                for j, seq in enumerate(batch_seqs):
                    seq_len = min(len(seq), max_len - 2)  # 考虑特殊token
                    
                    # 获取实际序列的嵌入（排除<cls>和<eos>）
                    seq_embedding = token_embeddings[j, 1:seq_len+1, :].mean(0)
                    all_embeddings.append(seq_embedding.cpu().numpy())
                    
                    # 简化的重要性评估（基于嵌入方差）
                    token_variance = token_embeddings[j, 1:seq_len+1, :].var(dim=1)
                    importance = token_variance.cpu().numpy()
                    all_attention_weights.append(importance)
                
                # 清理GPU内存
                del token_embeddings, results
                torch.cuda.empty_cache()
                
        return np.array(all_embeddings), all_attention_weights
    
    def get_embeddings_single(self, sequence: str) -> np.ndarray:
        """获取单个序列的嵌入（用于推理）"""
        data = [("seq", sequence)]
        batch_labels, batch_strs, batch_tokens = self.batch_converter(data)
        batch_tokens = batch_tokens.to(self.device)
        
        with torch.no_grad():
            results = self.model(batch_tokens, repr_layers=[36], return_contacts=False)
            token_embeddings = results["representations"][36]
            seq_len = len(sequence)
            seq_embedding = token_embeddings[0, 1:seq_len+1, :].mean(0)
            
        return seq_embedding.cpu().numpy()

class FunctionPredictor:
    """改进的功能预测器 - 优化版"""
    
    def __init__(self, esm_analyzer: ESMPeptideAnalyzer):
        self.esm_analyzer = esm_analyzer
        self.classifiers = {}
        self.feature_importances = {}
        self.esm_cache = {}  # 缓存ESM嵌入以加速
        
    def prepare_features_batch(self, sequences: List[str]) -> np.ndarray:
        """批量准备特征：ESM嵌入 + 理化性质 + 序列模式"""
        
        print(f"Preparing features for {len(sequences)} sequences...")
        
        # 1. ESM嵌入（批处理）
        esm_features, _ = self.esm_analyzer.get_embeddings_and_attention_batch(sequences)
        
        # 2. 理化性质特征
        print("Calculating physicochemical features...")
        phys_features = []
        for seq in tqdm(sequences, desc="Physicochemical"):
            phys_features.append(self._calculate_physicochemical(seq))
        phys_features = np.array(phys_features)
        
        # 3. Motif特征
        print("Calculating motif features...")
        motif_features = []
        for seq in tqdm(sequences, desc="Motifs"):
            motif_features.append(self._calculate_motif_features(seq))
        motif_features = np.array(motif_features)
        
        # 合并所有特征
        return np.hstack([esm_features, phys_features, motif_features])
    
    def prepare_features_single(self, sequence: str) -> np.ndarray:
        """为单个序列准备特征（用于推理）"""
        
        # 检查缓存
        if sequence in self.esm_cache:
            esm_features = self.esm_cache[sequence]
        else:
            esm_features = self.esm_analyzer.get_embeddings_single(sequence)
            # 限制缓存大小
            if len(self.esm_cache) < 1000:
                self.esm_cache[sequence] = esm_features
        
        phys_features = self._calculate_physicochemical(sequence)
        motif_features = self._calculate_motif_features(sequence)
        
        return np.hstack([esm_features, phys_features, motif_features])
    
    def _calculate_physicochemical(self, sequence: str) -> np.ndarray:
        """计算理化性质"""
        # 氨基酸性质字典
        hydrophobic = set('AILMFVPGW')
        positive = set('RKH')
        negative = set('DE')
        aromatic = set('FWY')
        polar = set('STNQCY')
        small = set('AGSVCT')
        
        length = len(sequence)
        if length == 0:
            return np.zeros(15)
            
        features = [
            length,
            sum(1 for aa in sequence if aa in hydrophobic) / length,
            sum(1 for aa in sequence if aa in positive) / length,
            sum(1 for aa in sequence if aa in negative) / length,
            sum(1 for aa in sequence if aa in aromatic) / length,
            sum(1 for aa in sequence if aa in polar) / length,
            sum(1 for aa in sequence if aa in small) / length,
            sequence.count('P') / length,  # Proline
            sequence.count('C') / length,  # Cysteine
            sequence.count('G') / length,  # Glycine
            # 净电荷
            (sum(1 for aa in sequence if aa in positive) - 
             sum(1 for aa in sequence if aa in negative)) / length,
            # 疏水矩
            self._hydrophobic_moment(sequence),
            # Boman指数（蛋白结合潜力）
            self._boman_index(sequence),
            # 不稳定性指数
            self._instability_index(sequence),
            # 等电点估计
            self._estimate_pi(sequence)
        ]
        
        return np.array(features)
    
    def _calculate_motif_features(self, sequence: str) -> np.ndarray:
        """计算功能motif特征"""
        motifs = {
            # AMP motifs
            'KLAK': 1 if 'KLAK' in sequence else 0,
            'RWR': 1 if 'RWR' in sequence else 0,
            'RRWW': 1 if 'RRWW' in sequence else 0,
            # AOP motifs  
            'HH': 1 if 'HH' in sequence else 0,
            'YY': 1 if 'YY' in sequence else 0,
            'HY': 1 if 'HY' in sequence else 0,
            # CPP motifs
            'RRR': 1 if 'RRR' in sequence else 0,
            'RKKR': 1 if 'RKKR' in sequence else 0,
            # AHP motifs
            'IPP': 1 if 'IPP' in sequence else 0,
            'VPP': 1 if 'VPP' in sequence else 0,
            'LPP': 1 if 'LPP' in sequence else 0,
            # 通用模式
            'KK': sequence.count('KK'),
            'RR': sequence.count('RR'),
            'WW': sequence.count('WW'),
            'PP': sequence.count('PP')
        }
        return np.array(list(motifs.values()))
    
    def _hydrophobic_moment(self, sequence: str) -> float:
        """计算疏水矩（α-螺旋）"""
        hydrophobicity = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }
        
        if len(sequence) < 3:
            return 0.0
            
        moment = 0.0
        angle = 100  # α-螺旋的角度
        
        for i, aa in enumerate(sequence[:min(11, len(sequence))]):
            if aa in hydrophobicity:
                h = hydrophobicity[aa]
                moment += h * np.cos(np.radians(i * angle))
                
        return abs(moment) / min(11, len(sequence))
    
    def _boman_index(self, sequence: str) -> float:
        """Boman指数：蛋白结合潜力"""
        boman_values = {
            'A': -0.53, 'R': 3.71, 'N': 1.91, 'D': 3.49, 'C': -3.31,
            'Q': 2.18, 'E': 2.84, 'G': 0.00, 'H': 1.44, 'I': -3.77,
            'L': -3.37, 'K': 4.11, 'M': -2.49, 'F': -4.22, 'P': -0.75,
            'S': 0.52, 'T': 0.07, 'W': -3.77, 'Y': -1.39, 'V': -3.26
        }
        
        if not sequence:
            return 0.0
            
        total = sum(boman_values.get(aa, 0) for aa in sequence)
        return total / len(sequence)
    
    def _instability_index(self, sequence: str) -> float:
        """简化的不稳定性指数"""
        unstable_pairs = ['CC', 'CW', 'WC', 'WW', 'MM', 'YY']
        score = sum(1 for pair in unstable_pairs if pair in sequence)
        return score / max(1, len(sequence) - 1)
    
    def _estimate_pi(self, sequence: str) -> float:
        """估计等电点"""
        positive = sequence.count('R') + sequence.count('K') + 0.5 * sequence.count('H')
        negative = sequence.count('D') + sequence.count('E')
        
        if not sequence:
            return 7.0
            
        charge_ratio = (positive - negative) / len(sequence)
        pi = 7.0 + charge_ratio * 3.5
        return max(2.0, min(14.0, pi))
    
    def train_improved_classifiers(self, sequences: List[str], labels: np.ndarray):
        """训练改进的分类器 - 优化版"""
        
        print("\nExtracting features with batch processing...")
        features = self.prepare_features_batch(sequences)
        
        # 保存特征用于后续分析
        np.save('extracted_features.npy', features)
        print("Saved features to extracted_features.npy")
        
        function_names = ['AMP', 'AOP', 'CPP', 'AHP']
        
        for i, func_name in enumerate(function_names):
            print(f"\n{'='*60}")
            print(f"Training {func_name} classifier...")
            print('='*60)
            
            func_labels = labels[:, i]
            
            # 统计正负样本
            pos_count = func_labels.sum()
            neg_count = len(func_labels) - pos_count
            print(f"Positive samples: {pos_count}")
            print(f"Negative samples: {neg_count}")
            
            if pos_count < 10:
                print(f"Warning: Too few positive samples for {func_name}, skipping...")
                continue
            
            # 分割数据 - 确保测试集有正样本
            X_train, X_test, y_train, y_test = train_test_split(
                features, func_labels, 
                test_size=0.2, 
                random_state=42, 
                stratify=func_labels if pos_count > 20 else None
            )
            
            # 训练随机森林 - 针对不平衡数据优化
            clf = RandomForestClassifier(
                n_estimators=300,  # 增加树的数量
                max_depth=30,
                min_samples_split=5,
                min_samples_leaf=2,
                class_weight='balanced_subsample',  # 更好的类权重策略
                random_state=42,
                n_jobs=-1,
                verbose=0
            )
            
            # 交叉验证
            if pos_count > 20:
                scores = cross_val_score(clf, X_train, y_train, cv=5, scoring='f1')
                print(f"Cross-validation F1: {scores.mean():.3f} (+/- {scores.std() * 2:.3f})")
            
            # 训练最终模型
            clf.fit(X_train, y_train)
            
            # 测试集评估
            y_pred = clf.predict(X_test)
            y_pred_proba = clf.predict_proba(X_test)[:, 1]
            
            # 计算多个指标
            from sklearn.metrics import roc_auc_score, precision_score, recall_score
            
            test_f1 = f1_score(y_test, y_pred)
            test_precision = precision_score(y_test, y_pred, zero_division=0)
            test_recall = recall_score(y_test, y_pred, zero_division=0)
            
            print(f"\nTest set performance:")
            print(f"  F1 Score: {test_f1:.3f}")
            print(f"  Precision: {test_precision:.3f}")
            print(f"  Recall: {test_recall:.3f}")
            
            if y_test.sum() > 0:  # 有正样本才计算AUC
                test_auc = roc_auc_score(y_test, y_pred_proba)
                print(f"  AUC: {test_auc:.3f}")
            
            # 保存分类器和特征重要性
            self.classifiers[func_name] = clf
            self.feature_importances[func_name] = clf.feature_importances_
            
            # 详细报告
            if func_name == 'AOP':  # 特别关注AOP
                print(f"\n{func_name} Detailed Report:")
                print(classification_report(y_test, y_pred))

class SequenceEditor:
    """序列编辑器：在保持原有功能的基础上添加新功能 - 优化版"""
    
    def __init__(self, esm_analyzer: ESMPeptideAnalyzer, predictor: FunctionPredictor):
        self.esm_analyzer = esm_analyzer
        self.predictor = predictor
        
    def identify_editable_positions(self, sequence: str, original_function: str) -> List[int]:
        """识别可编辑位置（不会破坏原功能的位置）"""
        
        # 简化版：基于氨基酸类型和功能需求
        protected_positions = set()
        
        if original_function == 'AMP':
            # 保护正电荷和关键疏水残基
            for i, aa in enumerate(sequence):
                if aa in 'RK' and i % 3 == 0:  # 每3个位置的RK
                    protected_positions.add(i)
                if aa in 'WFL' and i % 4 == 2:  # 疏水残基的特定模式
                    protected_positions.add(i)
                    
        elif original_function == 'CPP':
            # 保护所有精氨酸
            for i, aa in enumerate(sequence):
                if aa == 'R':
                    protected_positions.add(i)
                    
        elif original_function == 'AOP':
            # 保护抗氧化关键残基
            for i, aa in enumerate(sequence):
                if aa in 'HYWC':
                    protected_positions.add(i)
                    
        elif original_function == 'AHP':
            # 保护Pro位点
            for i in range(len(sequence)):
                if sequence[i] == 'P':
                    protected_positions.add(i)
                    if i > 0:
                        protected_positions.add(i-1)
        
        # 返回可编辑位置
        editable = [i for i in range(len(sequence)) if i not in protected_positions]
        
        # 至少保留30%可编辑
        min_editable = max(3, int(len(sequence) * 0.3))
        if len(editable) < min_editable:
            # 放松限制，只保护最关键的位置
            importance_scores = []
            for i in range(len(sequence)):
                if i in protected_positions:
                    if sequence[i] in 'RK' and original_function == 'AMP':
                        importance_scores.append((i, 10))
                    elif sequence[i] == 'R' and original_function == 'CPP':
                        importance_scores.append((i, 10))
                    else:
                        importance_scores.append((i, 5))
                        
            importance_scores.sort(key=lambda x: x[1], reverse=True)
            
            # 只保护最重要的位置
            protected_positions = set([pos for pos, score in importance_scores[:len(sequence)-min_editable]])
            editable = [i for i in range(len(sequence)) if i not in protected_positions]
            
        return editable
    
    def suggest_mutations(self, sequence: str, target_function: str, n_suggestions: int = 5) -> List[Dict]:
        """建议突变方案 - 优化版"""
        
        suggestions = []
        
        # 目标功能的优选氨基酸（按优先级排序）
        target_aas = {
            'AMP': ['K', 'R', 'W', 'L', 'F'],
            'AOP': ['H', 'Y', 'W', 'C', 'M'],
            'CPP': ['R', 'K', 'W'],
            'AHP': ['P', 'I', 'V', 'L', 'A']
        }
        
        # 获取当前序列的功能概率
        features = self.predictor.prepare_features_single(sequence)
        current_probs = {}
        for func, clf in self.predictor.classifiers.items():
            current_probs[func] = clf.predict_proba(features.reshape(1, -1))[0, 1]
        
        # 确定原始功能
        original_func = max(current_probs, key=current_probs.get)
        print(f"Original function: {original_func} (prob: {current_probs[original_func]:.3f})")
        
        # 识别可编辑位置
        editable_positions = self.identify_editable_positions(sequence, original_func)
        
        if not editable_positions:
            print(f"Warning: No editable positions found in {sequence}")
            return []
        
        print(f"Editable positions: {len(editable_positions)} / {len(sequence)}")
        
        # 策略1：单点突变
        for pos in editable_positions[:10]:  # 限制搜索空间
            for aa in target_aas[target_function][:3]:  # 只用最好的3个氨基酸
                if sequence[pos] == aa:
                    continue
                    
                mutant = list(sequence)
                mutant[pos] = aa
                mutant_seq = ''.join(mutant)
                
                # 评估突变体
                mut_features = self.predictor.prepare_features_single(mutant_seq)
                mut_probs = {}
                for func, clf in self.predictor.classifiers.items():
                    mut_probs[func] = clf.predict_proba(mut_features.reshape(1, -1))[0, 1]
                
                improvement = mut_probs[target_function] - current_probs[target_function]
                preservation = mut_probs[original_func] / max(0.01, current_probs[original_func])
                
                if improvement > 0.05:  # 降低阈值
                    suggestions.append({
                        'sequence': mutant_seq,
                        'mutations': [(pos, sequence[pos], aa)],
                        'improvement': improvement,
                        'preservation': preservation,
                        'new_probs': mut_probs,
                        'score': improvement * min(1.0, preservation)
                    })
        
        # 策略2：双点突变（如果单点不够）
        if len(suggestions) < n_suggestions and len(editable_positions) >= 2:
            import itertools
            
            for pos1, pos2 in itertools.combinations(editable_positions[:8], 2):
                for aa1, aa2 in itertools.product(target_aas[target_function][:2], repeat=2):
                    mutant = list(sequence)
                    mutant[pos1] = aa1
                    mutant[pos2] = aa2
                    mutant_seq = ''.join(mutant)
                    
                    mut_features = self.predictor.prepare_features_single(mutant_seq)
                    mut_probs = {}
                    for func, clf in self.predictor.classifiers.items():
                        mut_probs[func] = clf.predict_proba(mut_features.reshape(1, -1))[0, 1]
                    
                    improvement = mut_probs[target_function] - current_probs[target_function]
                    preservation = mut_probs[original_func] / max(0.01, current_probs[original_func])
                    
                    if improvement > 0.1 and preservation > 0.7:
                        suggestions.append({
                            'sequence': mutant_seq,
                            'mutations': [(pos1, sequence[pos1], aa1), (pos2, sequence[pos2], aa2)],
                            'improvement': improvement,
                            'preservation': preservation,
                            'new_probs': mut_probs,
                            'score': improvement * min(1.0, preservation)
                        })
                    
                    if len(suggestions) >= n_suggestions * 5:
                        break
        
        # 排序并返回最佳建议
        suggestions.sort(key=lambda x: x['score'], reverse=True)
        return suggestions[:n_suggestions]

def main():
    """主程序 - 优化版"""
    
    print("="*80)
    print("多功能肽智能改造系统 - 优化版")
    print("基于ESM-2-3B模型")
    print("="*80)
    
    # 1. 加载数据
    print("\n[1] Loading peptide data...")
    df = pd.read_csv('data/processed/training_data.csv')
    
    # 筛选长度合适的肽
    df_filtered = df[(df['length'] >= 10) & (df['length'] <= 30)].copy()
    
    sequences = df_filtered['sequence'].tolist()
    labels = df_filtered[['AMP', 'AOP', 'CPP', 'AHP']].values
    
    print(f"Loaded {len(sequences)} peptides (10-30 aa)")
    print(f"Function distribution:")
    for func in ['AMP', 'AOP', 'CPP', 'AHP']:
        count = df_filtered[func].sum()
        print(f"  {func}: {count} ({count/len(df_filtered)*100:.1f}%)")
    
    # 2. 初始化模型（优化批大小）
    print("\n[2] Initializing ESM-2 model...")
    esm_analyzer = ESMPeptideAnalyzer(batch_size=128)  # 对于96G显存，16是安全的
    
    # 3. 训练改进的预测器
    print("\n[3] Training improved function predictors...")
    predictor = FunctionPredictor(esm_analyzer)
    
    # 采样训练 - 使用更多数据
    sample_size = min(8000, len(sequences))  # 使用更多数据
    print(f"Sampling {sample_size} peptides for training...")
    
    if len(sequences) > sample_size:
        # 分层采样，确保每个功能都有足够样本
        sample_indices = []
        for func_idx in range(4):
            func_indices = np.where(labels[:, func_idx] == 1)[0]
            if len(func_indices) > 0:
                n_samples = min(len(func_indices), sample_size // 4)
                selected = np.random.choice(func_indices, n_samples, replace=False)
                sample_indices.extend(selected)
        
        # 填充剩余的
        remaining = set(range(len(sequences))) - set(sample_indices)
        if len(sample_indices) < sample_size and remaining:
            additional = np.random.choice(list(remaining), 
                                        min(sample_size - len(sample_indices), len(remaining)), 
                                        replace=False)
            sample_indices.extend(additional)
        
        sample_indices = np.array(sample_indices)
        train_sequences = [sequences[i] for i in sample_indices]
        train_labels = labels[sample_indices]
    else:
        train_sequences = sequences
        train_labels = labels
    
    predictor.train_improved_classifiers(train_sequences, train_labels)
    
    # 保存模型
    with open('improved_predictors_optimized.pkl', 'wb') as f:
        pickle.dump(predictor.classifiers, f)
    print("\nSaved improved predictors to improved_predictors_optimized.pkl")
    
    # 4. 序列改造实验
    print("\n[4] Sequence editing experiments...")
    editor = SequenceEditor(esm_analyzer, predictor)
    
    # 选择高质量的单功能肽进行改造
    single_func_peptides = df_filtered[df_filtered['num_functions'] == 1].copy()
    
    # 为每种功能选择测试案例
    test_cases = []
    for func in ['AMP', 'AOP', 'CPP', 'AHP']:
        func_peptides = single_func_peptides[single_func_peptides[func] == 1]
        if len(func_peptides) > 0:
            # 选择中等长度的肽
            candidates = func_peptides[(func_peptides['length'] >= 15) & 
                                      (func_peptides['length'] <= 25)]
            if len(candidates) > 0:
                selected = candidates.sample(min(5, len(candidates)))  # 每个功能选5个
                for _, row in selected.iterrows():
                    test_cases.append({
                        'sequence': row['sequence'],
                        'original_function': func,
                        'length': row['length']
                    })
    
    print(f"\nSelected {len(test_cases)} peptides for editing experiments")
    
    # 执行改造实验
    results = []
    successful_edits = 0
    
    for idx, case in enumerate(test_cases[:20]):  # 测试前20个
        seq = case['sequence']
        orig_func = case['original_function']
        
        print(f"\n{'='*60}")
        print(f"[{idx+1}/{min(20, len(test_cases))}] Original: {seq}")
        print(f"Function: {orig_func}, Length: {case['length']}")
        
        # 尝试添加每种其他功能
        for target_func in ['AMP', 'AOP', 'CPP', 'AHP']:
            if target_func == orig_func:
                continue
            
            print(f"\nTrying to add {target_func}...")
            suggestions = editor.suggest_mutations(seq, target_func, n_suggestions=3)
            
            if suggestions:
                best = suggestions[0]
                if best['score'] > 0.1:  # 只记录有意义的改进
                    print(f"  ✓ Success! Best mutation: {best['mutations']}")
                    print(f"    Improvement: {best['improvement']:.3f}")
                    print(f"    Preservation: {best['preservation']:.3f}")
                    print(f"    Score: {best['score']:.3f}")
                    
                    results.append({
                        'original_seq': seq,
                        'original_func': orig_func,
                        'target_func': target_func,
                        'mutant_seq': best['sequence'],
                        'mutations': str(best['mutations']),
                        'improvement': best['improvement'],
                        'preservation': best['preservation'],
                        'score': best['score'],
                        **{f'{k}_prob': v for k, v in best['new_probs'].items()}
                    })
                    successful_edits += 1
                else:
                    print(f"  ✗ No significant improvement found")
            else:
                print(f"  ✗ No valid mutations found")
    
    # 5. 保存和分析结果
    if results:
        results_df = pd.DataFrame(results)
        results_df.to_csv('editing_results_optimized.csv', index=False)
        print(f"\n[5] Results Summary")
        print("="*80)
        print(f"Total successful edits: {successful_edits}")
        print(f"Success rate: {successful_edits / (len(test_cases[:20]) * 3):.1%}")
        
        # 分析最成功的案例
        best_cases = results_df.nlargest(10, 'score')
        print("\nTop 10 successful edits:")
        print("-"*80)
        
        for idx, row in best_cases.iterrows():
            print(f"\n{idx+1}. {row['original_func']} → {row['original_func']}+{row['target_func']}")
            print(f"   Original: {row['original_seq'][:40]}...")
            print(f"   Mutant:   {row['mutant_seq'][:40]}...")
            print(f"   Mutations: {row['mutations']}")
            print(f"   Score: {row['score']:.3f} (Improvement: {row['improvement']:.3f}, Preservation: {row['preservation']:.3f})")
            print(f"   New probabilities:")
            for func in ['AMP', 'AOP', 'CPP', 'AHP']:
                if f'{func}_prob' in row:
                    print(f"     {func}: {row[f'{func}_prob']:.3f}")
        
        # 统计分析
        print("\n" + "="*80)
        print("Statistical Analysis:")
        print("-"*80)
        
        # 按功能组合分析成功率
        for orig in ['AMP', 'AOP', 'CPP', 'AHP']:
            orig_results = results_df[results_df['original_func'] == orig]
            if len(orig_results) > 0:
                print(f"\n{orig} as base function:")
                for target in ['AMP', 'AOP', 'CPP', 'AHP']:
                    if target != orig:
                        target_results = orig_results[orig_results['target_func'] == target]
                        if len(target_results) > 0:
                            avg_score = target_results['score'].mean()
                            max_score = target_results['score'].max()
                            print(f"  → {target}: Avg score: {avg_score:.3f}, Max: {max_score:.3f}")
    
    # 6. 生成可视化
    print("\n[6] Generating visualizations...")
    
    if results:
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. Score分布
        ax = axes[0, 0]
        results_df['score'].hist(bins=20, ax=ax, edgecolor='black')
        ax.set_xlabel('Edit Score')
        ax.set_ylabel('Count')
        ax.set_title('Distribution of Edit Scores')
        ax.axvline(0.5, color='red', linestyle='--', label='Success threshold')
        ax.legend()
        
        # 2. 功能改进热图
        ax = axes[0, 1]
        pivot_table = results_df.pivot_table(
            values='score',
            index='original_func',
            columns='target_func',
            aggfunc='mean'
        )
        sns.heatmap(pivot_table, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax)
        ax.set_title('Average Edit Score by Function Pair')
        
        # 3. Improvement vs Preservation散点图
        ax = axes[1, 0]
        scatter = ax.scatter(results_df['improvement'], 
                           results_df['preservation'],
                           c=results_df['score'], 
                           cmap='viridis',
                           s=50,
                           alpha=0.6,
                           edgecolors='black',
                           linewidth=0.5)
        ax.set_xlabel('Function Improvement')
        ax.set_ylabel('Original Function Preservation')
        ax.set_title('Trade-off between Improvement and Preservation')
        plt.colorbar(scatter, ax=ax, label='Score')
        
        # 4. 成功率条形图
        ax = axes[1, 1]
        success_rates = []
        labels_bar = []
        for orig in ['AMP', 'AOP', 'CPP', 'AHP']:
            orig_results = results_df[results_df['original_func'] == orig]
            if len(orig_results) > 0:
                success_rate = len(orig_results[orig_results['score'] > 0.3]) / max(1, len(orig_results))
                success_rates.append(success_rate)
                labels_bar.append(orig)
        
        if success_rates:
            ax.bar(labels_bar, success_rates)
            ax.set_ylabel('Success Rate')
            ax.set_xlabel('Original Function')
            ax.set_title('Edit Success Rate by Original Function')
            ax.set_ylim(0, 1)
            
            # 添加数值标签
            for i, (label, rate) in enumerate(zip(labels_bar, success_rates)):
                ax.text(i, rate + 0.02, f'{rate:.1%}', ha='center')
        
        plt.tight_layout()
        plt.savefig('editing_analysis_optimized.png', dpi=150, bbox_inches='tight')
        print("Saved analysis plots to editing_analysis_optimized.png")
    
    # 清理GPU内存
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        print(f"\nFinal GPU memory: {torch.cuda.memory_allocated() / 1e9:.2f} GB")
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)

if __name__ == "__main__":
    main()