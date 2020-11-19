# DetSimOptions
探测器模拟的主文件
### tut_detsim.py更新的代码
* 增加新光学模型的开关
```python
grp_pmt_op.add_argument("--new-optical-model", dest="new_optical_model", action="store_true",
                  help=mh("Use the new optical model."))
grp_pmt_op.add_argument("--old-optical-model", dest="new_optical_model", action="store_false",
                  help=mh("Use the old optical model"))
grp_pmt_op.set_defaults(new_optical_model=False)
```
* 
