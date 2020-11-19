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
* 如果想要找到与之相关的代码，请查询关键字new_optical_model
```python
if args.new_optical_model:
            op_process.property("UseAbsReemit").set(True)
            op_process.property("UseScintSimple").set(True)
        else:
            op_process.property("UseAbsReemit").set(False)
            op_process.property("UseScintSimple").set(False)
```

以及下面的代码
```python
 if args.new_optical_model:
            detsimfactory.property("GdLSAbsLengthMode").set(1)
        else:
            detsimfactory.property("GdLSAbsLengthMode").set(0)
```
